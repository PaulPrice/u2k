from collections import defaultdict
from logging import Logger

import numpy as np

from lsst.afw.geom import SkyWcs
from lsst.afw.image import ExposureF, ImageF, PhotoCalib, makeExposure, makeMaskedImage, FilterLabel
from lsst.afw.math import Warper
from lsst.daf.butler import Butler
from lsst.geom import Point2D
from lsst.meas.algorithms import CoaddPsf, CoaddPsfConfig
from lsst.pipe.tasks.coaddInputRecorder import CoaddInputRecorderTask
from lsst.skymap import BaseSkyMap, PatchInfo

from .detector import makeDetector
from .math import renormalizeVariance


def getWarpMapping(
    wcsList: list[SkyWcs],
    dimsList: list[tuple[int, int]],
    skymap: BaseSkyMap,
    sample: float = 500,
    log: Logger | None = None,
) -> dict[tuple[int, int], list[int]]:
    """Determine the mapping of patch to exposures
    
    Parameters
    ----------
    wcsList : `list[lsst.afw.geom.SkyWcs]`
        The WCSs for the exposures to warp.
    dimsList : `list[tuple[int, int]]`
        The dimensions of the exposures to warp (W, H).
    skymap : `lsst.skymap.SkyMap`
        The skymap to warp into.
    sample : `float`, optional
        The spacing at which to sample the input image when determining the
        overlap.
    log : `logging.Logger`, optional
        Logger to use for logging messages. If not provided, no messages will be
        logged.

    Returns
    -------
    patchToIndices : `dict[tuple[int, int], list[int]]`
        A mapping from patch (identified by tract and patch indices) to the list
        of indices of exposures that overlap with that patch.
    """
    if len(wcsList) != len(dimsList):
        raise ValueError("wcsList and dimsList must have the same length")
    exposureToPatches = [set() for _ in wcsList]
    for ii, (wcs, dims) in enumerate(zip(wcsList, dimsList)):
        width, height = dims
        xNumSamples = int(width / sample)
        yNumSamples = int(height / sample)
        xSample = np.linspace(0, width - 1, xNumSamples)
        ySample = np.linspace(0, height - 1, yNumSamples)
        coordList = [wcs.pixelToSky(Point2D(x, y)) for y in ySample for x in xSample]
        for tract, patches in skymap.findTractPatchList(coordList):
            for pp in patches:
                exposureToPatches[ii].add((tract.tract_id, pp.sequential_index))

    allPatches = set()
    for ii, patches in enumerate(exposureToPatches):
        allPatches.update(patches)
    if log is not None:
        log.info("Found %d patches across %d exposures", len(allPatches), len(wcsList))

    patchToIndices = defaultdict(list)
    for ii, patches in enumerate(exposureToPatches):
        for patch in patches:
            patchToIndices[patch].append(ii)

    return patchToIndices


def warpExposures(
    log: Logger,
    exposureList: list[ExposureF],
    patchInfo: PatchInfo,
    warpConfig: Warper.ConfigClass,
    *,
    butler: Butler | None = None,
    datasetType: str | None = None,
    dataId: dict | None = None,
    oldCollection: str | None = None,
) -> ExposureF:
    """Warp a list of exposures into a single patch of our skymap system

    Parameters
    ----------
    log : `logging.Logger`
        Logger to use for logging messages.
    exposureList : `list[lsst.afw.image.ExposureF]`
        The calibrated images to warp.
    patchInfo : `lsst.skymap.PatchInfo`
        The patch to warp into.
    warpConfig : `lsst.afw.math.Warper.ConfigClass`
        Configuration for the warping.
    butler : `lsst.daf.butler.Butler`, optional
        The butler to use for retrieving the old coadd if it exists. If not
        provided, the old coadd will not be included in the output.
    datasetType : `str`, optional
        The dataset type name for the old coadd.
    dataId : `dict`, optional
        The data ID for the old coadd.
    oldCollection : `str`, optional
        The collection containing the old coadd.

    Returns
    -------
    coadd : `lsst.afw.image.ExposureF`
        The warped coadd for the patch.
    """ 
    task = CoaddInputRecorderTask(log=log, name="inputRecorder")
    inputRecorder = task.makeCoaddTempExpRecorder(0, len(exposureList))

    warper = Warper.fromConfig(warpConfig)
    sumImage = ImageF(patchInfo.getOuterBBox())
    sumImage[:] = 0.0
    sumWeight = np.zeros_like(sumImage.array)
    totalGood = 0
    for exposure in exposureList:
        exposure.setDetector(makeDetector())

        warped = warper.warpExposure(
            patchInfo.getWcs(),
            exposure,
            maxBBox=patchInfo.getOuterBBox(),
            destBBox=patchInfo.getOuterBBox(),
        )

        good = np.isfinite(warped.image.array) & (warped.variance.array > 0)
        weight = 1.0/warped.variance.array[good]
        sumImage.array[good] += warped.image.array[good]*weight
        sumWeight[good] += weight
        numGood = np.sum(good)
        totalGood += numGood
        inputRecorder.addCalExp(exposure, 0, numGood)

    # Combine with the u2k v2 coadd if it exists.
    if butler is not None:
        dataRef = butler.find_dataset(datasetType, dataId, collections=oldCollection)
        if dataRef:
            old = butler.get(dataRef)
            assert old.getDimensions() == sumImage.getDimensions()
            old.maskedImage = old.getPhotoCalib().calibrateImage(old.maskedImage)
            old.setDetector(makeDetector())
            old.setFilter(exposureList[0].getFilter())
            good = np.isfinite(old.image.array) & (old.variance.array > 0)
            weight = 1.0/old.variance.array[good]
            sumImage.array[good] += old.image.array[good]*weight
            sumWeight[good] += weight
            numGood = np.sum(good)
            inputRecorder.addCalExp(old, 0, numGood)

    good = np.isfinite(sumWeight) & (sumWeight > 0)
    if good.sum() == 0:
        log.warning("Patch %s has no good pixels", dataId)
        return None
    log.info("Patch %s has %d good pixels", dataId, np.sum(good))

    sumImage.array[good] /= sumWeight[good]
    coadd = makeExposure(makeMaskedImage(sumImage))

    coadd.image.array[~good] = 0.0
    coadd.mask.array[~good] |= coadd.mask.getPlaneBitMask("NO_DATA")
    coadd.variance.array[good] = 1.0/sumWeight[good]
    coadd.variance.array[~good] = 0.0
    renormalizeVariance(coadd.maskedImage)

    # Set the photometric calibration
    coadd.setPhotoCalib(PhotoCalib(1.0))
    coadd.metadata["BUNIT"] = "nJy"

    coadd.setWcs(patchInfo.getWcs())
    coadd.setFilter(FilterLabel.fromBand(exposureList[0].getFilter().bandLabel))

    # Set the PSF and aperture corrections
    inputRecorder.finish(coadd, totalGood)
    coadd.setPsf(CoaddPsf(
        inputRecorder.coaddInputs.ccds,
        patchInfo.getWcs(),
        CoaddPsfConfig().makeControl(),
    ))

    return coadd
