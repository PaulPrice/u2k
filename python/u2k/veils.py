from __future__ import annotations

from collections import defaultdict
from logging import Logger

import numpy as np
import astropy.io.fits

from lsst.afw.image import ExposureF, ImageF, makeMaskedImage, makeExposure, PhotoCalib, FilterLabel
from lsst.afw.fits import readMetadata
from lsst.afw.math import Warper
from lsst.afw.geom import SkyWcs
from lsst.daf.butler import Butler
from lsst.geom import Point2D
from lsst.meas.algorithms import CoaddPsf, CoaddPsfConfig
from lsst.pipe.tasks.coaddInputRecorder import CoaddInputRecorderTask
from lsst.skymap import PatchInfo

from .butler import registerDatasetType
from .calibrateImage import CoaddCalibrateImageTask
from .detector import makeDetector
from .logging import getLogger
from .math import renormalizeVariance
from .processPool import ProcessPool


def readVeils(filename: str) -> ExposureF:
    """Read a VEILS image.

    Parameters
    ----------
    filename : `str`
        The path to the image file.
    """
    exposureList = []
    with astropy.io.fits.open(filename) as fits:
        for ii in range(16):
            image = fits[ii + 1].data.astype(np.float32)
            header = readMetadata(filename, hdu=ii + 1)
            image = ImageF(image)
            exposure = makeExposure(makeMaskedImage(image), SkyWcs(header))
            exposure.setMetadata(header)
            exposureList.append(exposure)
    return exposureList


def readVeilsConf(filename: str) -> list[ImageF]:
    """Read a VEILS confidence image.

    Parameters
    ----------
    filename : `str`
        The path to the confidence image file.
    """
    confList = []
    with astropy.io.fits.open(filename) as fits:
        for ii in range(16):
            conf = fits[ii + 1].data.astype(np.float32)
            confList.append(ImageF(conf))
    return confList


def warpVeilsPatch(
    log: Logger,
    butler: Butler,
    oldCollection: str,
    datasetType: str,
    dataId: dict[str, str|int],
    exposureList: list[ExposureF],
    patchInfo: PatchInfo,
    warpConfig: Warper.ConfigClass,
    overwrite: bool = False,
):
    """Warp VEILS images into a single patch of our skymap system
    
    Parameters
    ----------
    log : `logging.Logger`
        Logger to use for logging messages.
    butler : `lsst.daf.butler.Butler`
        The butler to use.
    oldCollection : `str`
        The collection with old data to check for existing output in.
    datasetType : `str`
        The dataset type name for the warped image.
    dataId : `dict`
        The data ID for the warped image.
    exposureList : list of `lsst.afw.image.ExposureF`
        The list of exposures to warp and coadd.
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

    task = CoaddInputRecorderTask(log=log, name="inputRecorder")
    inputRecorder = task.makeCoaddTempExpRecorder(0, len(exposureList))

    warper = Warper.fromConfig(warpConfig)
    image = ImageF(patchInfo.getOuterBBox())
    image[:] = 0.0
    sumWeight = np.zeros_like(image.array)
    totalGood = 0
    for exposure in exposureList:
        warped = warper.warpExposure(
            patchInfo.getWcs(),
            exposure,
            exposure.getWcs(),
            maxBBox=patchInfo.getOuterBBox(),
            destBBox=patchInfo.getOuterBBox(),
        )
        good = np.isfinite(warped.image.array) & (warped.variance.array > 0)
        weight = 1.0/warped.variance.array[good]
        image.array[good] += warped.image.array[good]*weight
        sumWeight[good] += weight
        numGood = np.sum(good)
        totalGood += numGood
        inputRecorder.addCalExp(exposure, 0, numGood)

    # Combine with the u2k v2 coadd if it exists.
    dataRef = butler.find_dataset(datasetType, dataId, collections=oldCollection)
    if dataRef:
        old = butler.get(dataRef)
        assert old.getDimensions() == image.getDimensions()
        old.maskedImage = old.getPhotoCalib().calibrateImage(old.maskedImage)
        old.setDetector(makeDetector())
        old.setFilter(exposureList[0].getFilter())
        good = np.isfinite(old.image.array) & (old.variance.array > 0)
        weight = 1.0/old.variance.array[good]
        image.array[good] += old.image.array[good]*weight
        sumWeight[good] += weight
        numGood = np.sum(good)
        inputRecorder.addCalExp(old, 0, numGood)

    good = np.isfinite(sumWeight) & (sumWeight > 0)
    log.info("Patch %s has %d good pixels", dataId, np.sum(good))
    image.array[good] /= sumWeight[good]
    coadd = makeExposure(makeMaskedImage(image))

    coadd.mask.array[~good] |= coadd.mask.getPlaneBitMask("NO_DATA")
    coadd.variance.array[good] = 1.0/sumWeight[good]
    coadd.variance.array[~good] = np.nan
    renormalizeVariance(coadd.maskedImage)

    # Set the photometric calibration
    coadd.setPhotoCalib(PhotoCalib(1.0))
    coadd.metadata["BUNIT"] = "nJy"

    # Set the PSF and aperture corrections
    inputRecorder.finish(coadd, totalGood)
    coadd.setPsf(CoaddPsf(
        inputRecorder.coaddInputs.ccds,
        patchInfo.getWcs(),
        CoaddPsfConfig().makeControl(),
    ))

    coadd.setWcs(patchInfo.getWcs())
    coadd.setFilter(exposureList[0].getFilter())

    butler.put(coadd, datasetType, dataId=dataId)


def calibrateVeils(
    exposure: ExposureF,
    conf: ImageF,
    log: Logger,
    calibrateConfig: CoaddCalibrateImageTask.ConfigClass,
) -> ExposureF:
    """Calibrate a VEILS image

    Parameters
    ----------
    dataId : any
        An identifier, used for logging purposes.
    exposure : `lsst.afw.image.ExposureF`
        The image to calibrate.
    conf : `lsst.afw.image.ImageF`
        The confidence image (inverse variance as a function of position).
    log : `logging.Logger`
        Logger to use for logging messages.
    calibrateConfig : `u2k.CoaddCalibrateImageTask.ConfigClass`
        Configuration for the calibration.

    Returns
    -------
    calibrated : `lsst.afw.image.ExposureF`
        The calibrated image.
    """
    good = np.isfinite(conf.array) & (conf.array > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        variance = 1.0/conf.array
    exposure.variance.array[:] = variance
    exposure.mask.array[~good] |= exposure.mask.getPlaneBitMask("NO_DATA")
    renormalizeVariance(exposure.maskedImage)

    zp = exposure.metadata["PHOTZP"]
    exposure.setPhotoCalib(PhotoCalib(3631e9 * 10**(-0.4*zp)))

    task = CoaddCalibrateImageTask(config=calibrateConfig, log=log)
    task.run(exposure=exposure)
    return exposure


def warpVeils(
    butler: Butler,
    exposureList: list[ExposureF],
    confList: list[ImageF],
    band: str,
    *,
    skymapName: str = "hsc",
    datasetType: str = "deepCoadd",
    oldCollection: str = "coadds",
    sample: float = 500,
    warpConfig: Warper.ConfigClass = Warper.ConfigClass(),
    calibrateConfig: CoaddCalibrateImageTask.ConfigClass = CoaddCalibrateImageTask.ConfigClass(),
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
    oldCollection : `str`, optional
        The collection with old data to check for existing output in.
    datasetType : `str`, optional
        The dataset type name for the warped image.    
    sample : `float`, optional
        Number of pixels to sample when determining the overlap.
    warpConfig : `lsst.afw.math.Warper.ConfigClass`, optional
        Configuration for the warping.
    calibrateConfig : `u2k.CoaddCalibrateImageTask.ConfigClass`, optional
        Configuration for the calibration.
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

    for exposure in exposureList:
        exposure.setFilter(FilterLabel.fromBandPhysical(band, band))

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

    pool = ProcessPool(workers)
    calibrated = pool.map(
        calibrateVeils,
        exposureList,
        confList,
        log=log,
        calibrateConfig=calibrateConfig,
    )

    with pool:
        for tractPatch in patchToIndices:
            tractId, patchId = tractPatch
            indices = patchToIndices[tractPatch]
            dataId = dict(instrument="u2k", band=band, tract=tractId, patch=patchId, skymap=skymapName)
            patchInfo = skymap[tractId][patchId]
            pool.submit(
                warpVeilsPatch,
                log,
                butler,
                oldCollection,
                datasetType,
                dataId,
                [calibrated[ii] for ii in indices],
                patchInfo,
                warpConfig,
                overwrite,
            )
