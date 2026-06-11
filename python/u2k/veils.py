from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from logging import Logger

import numpy as np
import astropy.io.fits

from lsst.afw.image import ExposureF, ImageF, makeMaskedImage, makeExposure, PhotoCalib
from lsst.afw.fits import readMetadata
from lsst.afw.math import Warper
from lsst.afw.geom import SkyWcs
from lsst.daf.butler import Butler
from lsst.geom import Point2D
from lsst.meas.algorithms.subtractBackground import SubtractBackgroundTask
from lsst.skymap import PatchInfo

from .butler import registerDatasetType
from .logging import getLogger
from .math import robustRms


def readVeils(filename: str) -> ExposureF:
    """Read a VEILS image.
    """
    with astropy.io.fits.open(filename) as fits:
        image = fits[1].data.astype(np.float32)
    header = readMetadata(filename)
    image = ImageF(image)
    exposure = makeExposure(makeMaskedImage(image), SkyWcs(header))
    exposure.setMetadata(header)
    return exposure


def warpVeilsPatch(
    log: Logger,
    butler: Butler,
    datasetType: str,
    dataId: dict[str, str|int],
    exposureList: list[ExposureF],
    confList: list[ImageF],
    patchInfo: PatchInfo,
    warpConfig: Warper.ConfigClass,
    overwrite: bool = False,
):
    """Warp an ASIAA image into a single patch of our skymap system
    
    Parameters
    ----------
    log : `logging.Logger`
        Logger to use for logging messages.
    butler : `lsst.daf.butler.Butler`
        The butler to use.
    datasetType : `str`
        The dataset type name for the warped image.
    dataId : `dict`
        The data ID for the warped image.
    exposureList : list of `lsst.afw.image.ExposureF`
        The list of exposures to warp and coadd.
    confList : list of `lsst.afw.image.ImageF`
        The list of confidence maps to warp.
    patchInfo : `lsst.skymap.PatchInfo`
        The patch to warp into.
    warpConfig : `lsst.afw.math.Warper.ConfigClass`
        Configuration for the warping.
    overwrite : `bool`, optional
        Whether to overwrite the output if it already exists.
    """
    dataRef = butler.find_dataset(datasetType, dataId)
    if dataRef:
        if overwrite:
            log.info("Overwriting existing output for %s", dataId)
            butler.pruneDatasets([dataRef], disassociate=True, unstore=True, purge=True)
        else:
            log.info("Output for %s already exists and overwrite=False, skipping", dataId)
            return

    log.info("Warping %d exposures into patch %s", len(exposureList), dataId)

    warper = Warper.fromConfig(warpConfig)
    image = ImageF(patchInfo.getOuterBBox())
    image[:] = 0.0
    sumConf = np.zeros_like(image.array)
    for exposure, conf in zip(exposureList, confList):
        zp = exposure.metadata["PHOTZP"]
        photoCalib = PhotoCalib(3631e9 * 10**(-0.4*zp))
        calibrated = photoCalib.calibrateImage(exposure.maskedImage)

        warpedImage = warper.warpImage(
            patchInfo.getWcs(),
            calibrated.image,
            exposure.getWcs(),
            maxBBox=patchInfo.getOuterBBox(),
            destBBox=patchInfo.getOuterBBox(),
        )
        warpedConf = warper.warpImage(
            patchInfo.getWcs(),
            conf,
            exposure.getWcs(),
            maxBBox=patchInfo.getOuterBBox(),
            destBBox=patchInfo.getOuterBBox(),
        )
        good = np.isfinite(warpedImage.array) & (warpedConf.array > 0)
        image.array[good] += warpedImage.array[good]*warpedConf.array[good]
        sumConf[good] += warpedConf.array[good]

    good = np.isfinite(sumConf) & (sumConf > 0)
    log.info("Patch %s has %d good pixels", dataId, np.sum(good))
    image.array[good] /= sumConf[good]
    exposure = makeExposure(makeMaskedImage(image))

    exposure.mask.array[~good] |= exposure.mask.getPlaneBitMask("NO_DATA")
    exposure.variance.array[good] = 1.0/sumConf[good]
    exposure.variance.array[~good] = np.nan

    # Calibrate the variance plane
    sigNoise = exposure.image.array/np.sqrt(exposure.variance.array)
    rms = robustRms(sigNoise[good])
    exposure.variance.array *= rms**2

    # Set the photometric calibration
    exposure.setPhotoCalib(PhotoCalib(1.0))
    exposure.metadata["BUNIT"] = "nJy"
    butler.put(exposure, datasetType, dataId=dataId)


def warpVeils(
    butler: Butler,
    exposureList: list[ExposureF],
    confList: list[ImageF],
    band: str,
    *,
    skymapName: str = "hsc",
    datasetType: str = "deepCoadd_warp",
    sample: float = 500,
    warpConfig: Warper.ConfigClass = Warper.ConfigClass(),
    backgroundConfig: SubtractBackgroundTask.ConfigClass = SubtractBackgroundTask.ConfigClass(),
    workers: int = 1,
    overwrite: bool = False,
):
    """Warp VEILS images into our skymap system
    
    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        The butler to use.
    exposure : `lsst.afw.image.ExposureF`
        The image to warp.
    timeImage : `lsst.afw.image.ImageF`
        The time image (exposure time as a function of position). This informs
        the variance plane.
    band : `str`
        The band name for the image.
    skymapName : `str`, optional
        The skymap name defining the output WCSes.
    datasetType : `str`, optional
        The dataset type name for the warped image.    
    sample : `float`, optional
        Number of pixels to sample when determining the overlap.
    warpConfig : `lsst.afw.math.Warper.ConfigClass`, optional
        Configuration for the warping.
    backgroundConfig : `lsst.meas.algorithms.subtractBackground.SubtractBackgroundTask.ConfigClass`, optional
        Configuration for the background subtraction.
    workers : `int`, optional
        Number of worker processes to use for warping the patches.
    overwrite : `bool`, optional
        Whether to overwrite the output if it already exists.
    """
    if len(exposureList) != len(confList):
        raise ValueError("exposureList and confList must have the same length")
    log = getLogger("u2k.warpVeils")
    registerDatasetType(butler.registry, datasetType, ("band", "tract", "patch", "skymap"), "ExposureF")
    skymap = butler.get("skyMap", collections="skymaps", skymap=skymapName)

    # Background subtraction
    bgTask = SubtractBackgroundTask(config=backgroundConfig, log=log)
    for exposure, conf in zip(exposureList, confList):
        good = np.isfinite(conf.array) & (conf.array > 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            variance = 1.0/conf.array
        exposure.variance.array[:] = variance
        exposure.mask.array[~good] |= exposure.mask.getPlaneBitMask("NO_DATA")

        with np.errstate(divide="ignore", invalid="ignore"):
            sigNoise = exposure.image.array/np.sqrt(exposure.variance.array)
        rms = robustRms(sigNoise[good])
        exposure.variance.array *= rms**2

        bgTask.run(exposure)

    exposureToPatches = [set() for _ in exposureList]
    for ii, exposure in enumerate(exposureList):
        wcs = exposure.getWcs()
        width, height = exposure.getDimensions()
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
    log.info("Found %d patches across %d exposures", len(allPatches), len(exposureList))

    patchToIndices = defaultdict(list)
    for ii, patches in enumerate(exposureToPatches):
        for patch in patches:
            patchToIndices[patch].append(ii)

    executor = ProcessPoolExecutor(workers)
    jobs = []
    for tractPatch in patchToIndices:
        tractId, patchId = tractPatch
        indices = patchToIndices[tractPatch]

        dataId = dict(instrument="u2k", band=band, tract=tractId, patch=patchId, skymap=skymapName)
        patchInfo = skymap[tractId][patchId]
        jobs.append(
            executor.submit(
                warpVeilsPatch,
                log,
                butler,
                datasetType,
                dataId,
                [exposureList[ii] for ii in indices],
                [confList[ii] for ii in indices],
                patchInfo,
                warpConfig,
                overwrite,
            )
        )

    try:
        for job in jobs:
            job.result()
    except Exception as exc:
        log.error("Error warping patch: %s", exc)
        for jj in jobs:
            jj.cancel()
        raise
    finally:
        executor.shutdown(wait=True)
