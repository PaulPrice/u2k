import numpy as np
from scipy.stats import qmc

import lsst.afw.detection
import lsst.afw.geom
import lsst.meas.algorithms.skyObjects
import lsst.pipe.tasks.multiBand

from lsst.afw.image import PhotoCalib

__all__ = ("DetectCoaddSourcesTask", "generateSkyObjects")


class DetectCoaddSourcesTask(lsst.pipe.tasks.multiBand.DetectCoaddSourcesTask):
    def run(self, exposure, *args, **kwargs):
        """Strip ``DETECTED`` mask plane before running detection

        We're operating on ``deepCoadd_calexp`` exposures from the U2K version 2
        processing; these have already had detection run on them. The old
        ``DETECTED`` mask plane is getting in the way of our new detection, so
        we clear it first.

        We also need to rescale the image so that it's in units of nJy (not with
        a magnitude zero point of 27 as previous LSST versions used), as this
        is what the TransformObjectCatalogTask expects.
        """
        exposure.mask.array &= ~exposure.mask.getPlaneBitMask("DETECTED")

        photoCalib = exposure.getPhotoCalib()
        assert photoCalib is not None, "Exposure has no PhotoCalib"

        if photoCalib.getCalibrationMean() == 27.0:
            # A magnitude zero-point is getting confused as a flux zero-point
            # We want a ZP of 31.4 mag for nJy, and we have 27, so a mag difference of 4.4
            photoCalib = PhotoCalib(10**(0.4*4.4))

        exposure.setMaskedImage(photoCalib.calibrateImage(exposure.maskedImage))
        exposure.setPhotoCalib(PhotoCalib(1.0))

        return super().run(exposure, *args, **kwargs)


def generateSkyObjects(mask, seed, config):
    """Generate a list of Footprints of sky objects

    Sky objects don't overlap with other objects. This is determined
    through the provided `mask` (in which objects are typically flagged
    as `DETECTED`).

    Sky objects are positioned using a quasi-random Halton sequence number
    generator. This is a deterministic sequence that mimics a random trial and
    error approach whilst acting to minimize clustering of points for a given
    field of view. Up to `nTrialSources` points are generated, returning the
    first `nSources` that do not overlap with the mask.

    Parameters
    ----------
    mask : `lsst.afw.image.Mask`
        Input mask plane, which identifies pixels to avoid for the sky
        objects.
    seed : `int`
        Random number generator seed.
    config : `SkyObjectsConfig`
        Configuration for finding sky objects.

    Returns
    -------
    skyFootprints : `list` of `lsst.afw.detection.Footprint`
        Footprints of sky objects. Each will have a peak at the center
        of the sky object.
    """
    if config.nSources <= 0:
        return []

    skySourceRadius = config.sourceRadius
    nSkySources = config.nSources
    nTrialSkySources = config.nTrialSources
    if nTrialSkySources is None:
        nTrialSkySources = config.nTrialSourcesMultiplier*nSkySources

    box = mask.getBBox()
    box.grow(-(int(skySourceRadius) + 1))  # Avoid objects partially off the image
    xMin, yMin = box.getMin()
    xMax, yMax = box.getMax()

    # Find bounds of unmasked pixels
    maskVal = mask.getPlaneBitMask("NO_DATA")
    unmasked = (mask.getArray() & maskVal) == 0
    unmaskedColumns = np.any(unmasked, axis=0)
    unmaskedRows = np.any(unmasked, axis=1)
    x0, y0 = mask.getXY0()
    xMin = max(xMin, x0 + int(np.where(unmaskedColumns)[0][0]))
    xMax = min(xMax, x0 + int(np.where(unmaskedColumns)[0][-1]))
    yMin = max(yMin, y0 + int(np.where(unmaskedRows)[0][0]))
    yMax = min(yMax, y0 + int(np.where(unmaskedRows)[0][-1]))

    avoid = lsst.afw.geom.SpanSet.fromMask(mask, mask.getPlaneBitMask(config.avoidMask))
    if config.growMask > 0:
        avoid = avoid.dilated(config.growMask)

    sampler = qmc.Halton(d=2, seed=seed).random(nTrialSkySources)
    sample = qmc.scale(sampler, [xMin, yMin], [xMax, yMax])

    skyFootprints = []
    for x, y in zip(sample[:, 0].astype(int), sample[:, 1].astype(int)):
        if len(skyFootprints) == nSkySources:
            break

        spans = lsst.afw.geom.SpanSet.fromShape(int(skySourceRadius), offset=(x, y))
        if spans.overlaps(avoid):
            continue

        fp = lsst.afw.detection.Footprint(spans, mask.getBBox())
        fp.addPeak(x, y, 0)
        skyFootprints.append(fp)

        # Add doubled-in-size sky object spanSet to the avoid mask.
        avoid = avoid.union(spans.dilated(int(skySourceRadius)))

    return skyFootprints


# Monkey-patch our improved function into place
lsst.meas.algorithms.skyObjects.generateSkyObjects = generateSkyObjects
