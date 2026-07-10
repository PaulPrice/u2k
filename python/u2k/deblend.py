# This file is part of meas_extensions_scarlet.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

import logging
import time
import traceback
from typing import Any

import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.scarlet.lite as scl
import numpy as np

from lsst.meas.extensions.scarlet import utils
from lsst.meas.extensions.scarlet.scarletDeblendTask import (
    DeblenderError,
    DeblenderSkippedError,
    ScarletDeblendConfig,
    ScarletDeblendTask,
    deblend,
    isPseudoSource,
)

__all__ = [
    "PiecewiseScarletDeblendConfig",
    "PiecewiseScarletDeblendTask",
]

logger = logging.getLogger(__name__)


class _PiecewiseTiledBlend:
    """Combines sources from multiple tile deblends into one blend-like object.

    This presents the same interface as ``scl.Blend`` as consumed by
    `ScarletDeblendTask._deblendParent`: it exposes ``sources``, ``loss``,
    ``components``, ``get_model()``, and ``to_data()``.

    Parameters
    ----------
    sources :
        The owned scarlet sources collected across all tiles — one entry per
        peak whose primary-zone tile they belong to.
    loss :
        Concatenated per-iteration loss values from all tile fits.
    stitchedModel :
        The model image built by summing every owned source's model.
    footprintBbox :
        The scarlet ``Box`` of the full (un-tiled) footprint, used as the
        bounding box reported in ``to_data()``.
    """

    def __init__(
        self,
        sources: list,
        loss: list[float],
        stitchedModel: scl.Image,
        footprintBbox: scl.Box,
    ):
        self.sources = sources
        self.loss = loss
        self._stitchedModel = stitchedModel
        self._bbox = footprintBbox

    @property
    def bbox(self) -> scl.Box:
        return self._bbox

    @property
    def components(self) -> list:
        return [c for src in self.sources for c in src.components]

    def get_model(self, convolve: bool = False, use_flux: bool = False) -> scl.Image:
        return self._stitchedModel

    def to_data(self) -> scl.io.ScarletBlendData:
        """Serialise the piecewise result as a single ``ScarletBlendData``.

        The sources dict is keyed by their catalog record ID (stored in
        ``source.metadata["id"]`` after `_addDeblendedSource` is called),
        exactly mirroring the behaviour of ``scl.Blend.to_data()``.
        """
        sources: dict[Any, scl.io.ScarletSourceBaseData] = {}
        for sidx, source in enumerate(self.sources):
            metadata = source.metadata or {}
            if "id" in metadata:
                sources[metadata["id"]] = source.to_data()
            else:
                sources[sidx] = source.to_data()
        return scl.io.ScarletBlendData(
            origin=self._bbox.origin,  # type: ignore[arg-type]
            shape=self._bbox.shape,   # type: ignore[arg-type]
            sources=sources,
            metadata=None,
        )


class PiecewiseScarletDeblendConfig(ScarletDeblendConfig):
    """Configuration for `PiecewiseScarletDeblendTask`.

    Inherits all fields from `ScarletDeblendConfig` and adds two parameters
    that control the tiling of oversized footprints.
    """

    boxSize = pexConfig.Field[int](
        default=300,
        doc=(
            "Side length (pixels) of the primary region of each tile used when "
            "breaking up a large footprint for piecewise deblending. "
            "The actual deblending box for each tile is larger by ``bufferSize`` "
            "on each internal edge."
        ),
    )
    bufferSize = pexConfig.Field[int](
        default=50,
        doc=(
            "Width (pixels) of the overlap buffer added on each *internal* tile "
            "boundary. Sources whose peaks lie near a division boundary are "
            "therefore modelled with their neighbours present, reducing edge "
            "artefacts. External edges of the parent footprint carry no buffer."
        ),
    )

    def setDefaults(self):
        super().setDefaults()
        self.maxFootprintSize = self.boxSize


class PiecewiseScarletDeblendTask(ScarletDeblendTask):
    """Variant of `ScarletDeblendTask` that deblends large footprints
    piecewise.

    Instead of skipping footprints that exceed the size thresholds configured
    in `ScarletDeblendConfig` (``maxFootprintArea``, ``maxFootprintSize``,
    ``maxAreaTimesPeaks``), this task divides them into a 2-D grid of
    overlapping tiles, deblends each tile independently with scarlet, and
    stitches the results together as if the whole footprint had been deblended
    in one pass.

    Each tile has a *primary region* of size ``boxSize × boxSize`` and a
    *deblending box* that extends ``bufferSize`` pixels beyond each internal
    boundary.  A peak is *owned* by the tile whose primary region contains
    that peak's centre; buffer-zone peaks are modelled but their deblended
    sources are discarded.  The stitched model is the sum of all owned
    sources' models.

    Small footprints (those that would pass the monolithic size checks) are
    passed directly to the parent `ScarletDeblendTask` logic unchanged.
    """

    ConfigClass = PiecewiseScarletDeblendConfig
    _DefaultName = "piecewiseScarletDeblend"

    # ------------------------------------------------------------------
    # Size-check overrides
    # ------------------------------------------------------------------

    def _isLargeFootprint(self, footprint: afwDet.Footprint) -> bool:
        """Always return False: piecewise task never skips for size."""
        return False

    def _isTooLargeForMonolithic(self, footprint: afwDet.Footprint) -> bool:
        """Return True if the parent class would skip this footprint for
        size."""
        return ScarletDeblendTask._isLargeFootprint(self, footprint)

    # ------------------------------------------------------------------
    # Override of _deblendParent to route large footprints piecewise
    # ------------------------------------------------------------------

    def _deblendParent(
        self,
        blendRecord: afwTable.SourceRecord,
    ) -> tuple[Any, scl.Image, scl.Image]:
        footprint = blendRecord.getFootprint()
        if self._isTooLargeForMonolithic(footprint):
            return self._deblendParentPiecewise(blendRecord)
        return super()._deblendParent(blendRecord)

    # ------------------------------------------------------------------
    # Tiling helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _tileAxis(
        lo: int,
        hi: int,
        boxSize: int,
        bufferSize: int,
    ) -> list[tuple[int, int, int, int]]:
        """Divide the interval ``[lo, hi)`` into overlapping tiles.

        Parameters
        ----------
        lo, hi :
            Half-open interval to tile (pixel coordinates, exclusive ``hi``).
        boxSize :
            Width of each tile's primary region.
        bufferSize :
            Extra pixels added on each *internal* edge.

        Returns
        -------
        tiles :
            List of ``(primaryLo, primaryHi, boxLo, boxHi)`` tuples, where
            the primary interval is half-open ``[primaryLo, primaryHi)`` and
            the deblending box is half-open ``[boxLo, boxHi)``.
        """
        tiles = []
        pos = lo
        while pos < hi:
            primaryLo = pos
            primaryHi = min(pos + boxSize, hi)
            # No buffer on the external edges of the full footprint
            boxLo = primaryLo if primaryLo == lo else primaryLo - bufferSize
            boxHi = primaryHi if primaryHi == hi else primaryHi + bufferSize
            tiles.append((primaryLo, primaryHi, boxLo, boxHi))
            pos = primaryHi
        return tiles

    def _makeTileFootprint(
        self,
        parentFootprint: afwDet.Footprint,
        tileBox: geom.Box2I,
    ) -> afwDet.Footprint | None:
        """Build a sub-footprint for one tile.

        The spans are the intersection of the parent footprint spans with
        ``tileBox``.  All peaks from the parent footprint whose centres fall
        inside ``tileBox`` are copied across (including buffer-zone peaks, so
        that scarlet models them as neighbours even if we don't retain them).

        Parameters
        ----------
        parentFootprint :
            The full deconvolved parent footprint.
        tileBox :
            The deblending bounding box for this tile (inclusive on both ends
            in ``Box2I`` convention).

        Returns
        -------
        tileFoot :
            The tile footprint, or ``None`` if the intersection is empty.
        """
        # Build a mask image that covers the tile box completely (all 1s).
        tileArray = np.ones(
            (tileBox.getHeight(), tileBox.getWidth()), dtype=np.int32
        )
        tileMask = afwImage.MaskX(tileArray, xy0=tileBox.getMin())

        # Intersect the parent spans with the tile box.
        tileSpans = parentFootprint.spans.intersect(tileMask, 1)
        if tileSpans.getArea() == 0:
            return None

        tileFoot = afwDet.Footprint(tileSpans)
        # Preserve the full peak schema (including merge_peak_* fields) so
        # that peak records assigned via peakSchemaMapper later are compatible.
        tileFoot.setPeakSchema(parentFootprint.peaks.getSchema())

        # Copy every peak whose centre falls inside the tile box.
        for peak in parentFootprint.peaks:
            if tileBox.contains(geom.Point2I(peak.getIx(), peak.getIy())):
                tileFoot.peaks.append(peak)

        return tileFoot

    # ------------------------------------------------------------------
    # Piecewise deblending
    # ------------------------------------------------------------------

    def _deblendParentPiecewise(
        self,
        blendRecord: afwTable.SourceRecord,
    ) -> tuple[_PiecewiseTiledBlend, scl.Image, scl.Image]:
        """Deblend a large parent footprint by tiling.

        The footprint is divided into a 2-D grid of overlapping tiles.  Each
        tile is deblended independently; only the sources whose peaks fall in
        the tile's *primary region* (the un-buffered interior) are retained.
        The results are stitched into a combined model image and returned as a
        ``_PiecewiseTiledBlend``.

        Parameters
        ----------
        blendRecord :
            The deconvolved parent source record to deblend.

        Returns
        -------
        combinedBlend :
            A ``_PiecewiseTiledBlend`` containing all owned sources.
        stitchedModel :
            The model image (sum of owned sources' models).
        chi2 :
            Per-pixel chi² image computed on the stitched model.

        Raises
        ------
        DeblenderSkippedError
            If the footprint should be skipped (pseudo-source or masked).
        DeblenderError
            If all tiles fail or piecewise deblending otherwise cannot proceed.
        """
        footprint = blendRecord.getFootprint()
        peaks = footprint.getPeaks()

        # Propagate the first peak's flags to the blend record (same as
        # base class).
        blendRecord.assign(peaks[0], self.parentPeakSchemaMapper)

        # Apply skip checks that are still relevant (pseudo-source, masked).
        # We deliberately do NOT check _isLargeFootprint (overridden to False)
        # or maxNumberOfPeaks (applied per tile below instead).
        if isPseudoSource(blendRecord, self.config.pseudoColumns):
            raise DeblenderSkippedError(
                f"Parent {blendRecord.getId()} is a pseudo source",
                blendRecord.getId(),
                "deblend_skipped_isPseudo",
            )
        if self._isMasked(footprint, self.mExposure):
            raise DeblenderSkippedError(
                f"Parent {blendRecord.getId()}: skipping masked footprint",
                blendRecord.getId(),
                "deblend_skipped_masked",
            )

        config = self.config
        bbox = footprint.getBBox()
        x0, y0 = bbox.getMinX(), bbox.getMinY()
        # _tileAxis uses half-open intervals, so pass getMax() + 1.
        x1, y1 = bbox.getMaxX() + 1, bbox.getMaxY() + 1

        xTiles = self._tileAxis(x0, x1, config.boxSize, config.bufferSize)
        yTiles = self._tileAxis(y0, y1, config.boxSize, config.bufferSize)

        nBands = len(self.context.observation.bands)
        width, height = bbox.getDimensions()

        # Initialise the stitched model covering the full footprint extent.
        stitchedModel = scl.Image(
            np.zeros((nBands, height, width), dtype=self.mExposure.image.array.dtype),
            bands=self.context.observation.images.bands,
            yx0=(y0, x0),
        )

        allOwnedSources: list = []
        allLoss: list[float] = []

        t0 = time.monotonic()

        for yPrimLo, yPrimHi, yBoxLo, yBoxHi in yTiles:
            for xPrimLo, xPrimHi, xBoxLo, xBoxHi in xTiles:
                # Inclusive Box2I uses max = exclusive_hi - 1.
                tileBox = geom.Box2I(
                    geom.Point2I(xBoxLo, yBoxLo),
                    geom.Point2I(xBoxHi - 1, yBoxHi - 1),
                )
                primaryBox = geom.Box2I(
                    geom.Point2I(xPrimLo, yPrimLo),
                    geom.Point2I(xPrimHi - 1, yPrimHi - 1),
                )

                tileFoot = self._makeTileFootprint(footprint, tileBox)
                if tileFoot is None or len(tileFoot.peaks) == 0:
                    continue

                # Per-tile peak count guard (replaces the footprint-level check
                # we bypassed above).
                if (
                    config.maxNumberOfPeaks > 0
                    and len(tileFoot.peaks) > config.maxNumberOfPeaks
                ):
                    self.log.warning(
                        "Tile (%d, %d) of parent %d has %d peaks > maxNumberOfPeaks=%d; "
                        "skipping tile",
                        xBoxLo, yBoxLo, blendRecord.getId(),
                        len(tileFoot.peaks), config.maxNumberOfPeaks,
                    )
                    continue

                # Per-tile spectrum-initialisation decision (the tile is
                # smaller than the full footprint, so this is more often
                # True).
                if config.setSpectra:
                    if config.maxSpectrumCutoff <= 0:
                        tileSpectrumInit = True
                    else:
                        tileSpectrumInit = (
                            len(tileFoot.peaks) * tileBox.getArea()
                            < config.maxSpectrumCutoff
                        )
                else:
                    tileSpectrumInit = False

                try:
                    tileBlend = deblend(
                        self.context, tileFoot, config, tileSpectrumInit
                    )
                except Exception as e:
                    blendError = type(e).__name__
                    if config.catchFailures:
                        self.log.warning(
                            "Error deblending tile (%d, %d) of parent %d: %s",
                            xBoxLo, yBoxLo, blendRecord.getId(), e,
                        )
                        traceback.print_exc()
                        continue
                    raise DeblenderError(
                        f"Unable to deblend tile of parent {blendRecord.getId()}: "
                        f"{blendError}",
                        blendRecord.getId(),
                        blendError,
                    )

                allLoss.extend(tileBlend.loss)

                # Retain only sources whose peak is in this tile's primary
                # zone.
                for scarletSource in tileBlend.sources:
                    peak = scarletSource.detectedPeak
                    if primaryBox.contains(geom.Point2I(peak.getIx(), peak.getIy())):
                        allOwnedSources.append(scarletSource)
                        # Add this source's model contribution to the stitch.
                        stitchedModel.insert(scarletSource.get_model())

        if not allOwnedSources:
            raise DeblenderError(
                f"Piecewise deblending of parent {blendRecord.getId()} "
                f"yielded no sources",
                blendRecord.getId(),
                "NoSourcesError",
            )

        tf = time.monotonic()
        runtime = (tf - t0) * 1000

        # Compute chi2 on the stitched model (non-zero pixels of the model
        # define the footprint mask for this calculation).
        stitchedFootprintMask = stitchedModel.data > 0
        chi2 = utils.calcChi2(
            stitchedModel, self.context.observation, stitchedFootprintMask
        )

        nComponents = sum(len(s.components) for s in allOwnedSources)
        nChild = len(allOwnedSources)

        # Convergence: use the last two loss values across all tiles as a
        # proxy. convergenceFailed=True means the fit ran but did not converge.
        if len(allLoss) >= 2:
            convergenceFailed = not (
                abs(allLoss[-2] - allLoss[-1])
                < config.relativeError * abs(allLoss[-1])
            )
        else:
            # No fitting ran across any tile; treat as vacuously converged.
            convergenceFailed = False

        # Provide at least two loss values so callers don't crash on len
        # checks.
        loss = allLoss if len(allLoss) >= 2 else [0.0, 0.0]

        combinedBlend = _PiecewiseTiledBlend(
            sources=allOwnedSources,
            loss=loss,
            stitchedModel=stitchedModel,
            footprintBbox=utils.bboxToScarletBox(bbox),
        )

        nPixels = int(np.sum(stitchedFootprintMask))
        self._updateParentRecord(
            parentRecord=blendRecord,
            nPeaks=len(peaks),
            nChild=nChild,
            nComponents=nComponents,
            runtime=runtime,
            iterations=len(allLoss),
            logL=allLoss[-1] if allLoss else np.nan,
            chi2=float(np.sum(chi2.data)) / max(1, nPixels),
            spectrumInit=False,  # varies per tile; report False as conservative
            convergenceFailed=convergenceFailed,
        )

        return combinedBlend, stitchedModel, chi2
