from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from logging import Logger

import numpy as np

from lsst.afw.image import ExposureF, ImageF, makeMaskedImage, makeExposure, PhotoCalib
from lsst.afw.math import Warper
from lsst.daf.butler import Butler
from lsst.geom import Point2D, Point2I, Box2I, SpherePoint
from lsst.skymap import PatchInfo

from .butler import registerDatasetType
from .logging import getLogger
from .math import robustRms


def warpAsiaaPatch(
    log: Logger,
    butler: Butler,
    datasetType: str,
    dataId: dict[str, str|int],
    exposure: ExposureF,
    timeImage: ImageF,
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
    image : `lsst.afw.image.ExposureF`
        The image to warp.
    timeImage : `lsst.afw.image.ImageF`
        The time image (exposure time as a function of position). This informs
        the variance plane.
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

    log.info("Warping patch %s", dataId)

    warper = Warper.fromConfig(warpConfig)
    warpedImage = warper.warpImage(
        patchInfo.getWcs(),
        exposure.image,
        exposure.getWcs(),
        maxBBox=patchInfo.getOuterBBox(),
        destBBox=patchInfo.getOuterBBox(),
    )
    warpedExposure = makeExposure(makeMaskedImage(warpedImage))

    # The variance goes as the background, which goes as the exposure time
    warpedExposure.variance[:] = warper.warpImage(
        patchInfo.getWcs(),
        timeImage,
        exposure.getWcs(),
        maxBBox=patchInfo.getOuterBBox(),
        destBBox=patchInfo.getOuterBBox(),
    )

    good = np.isfinite(warpedExposure.image.array)
    good &= np.isfinite(warpedExposure.variance.array) & (warpedExposure.variance.array > 0)

    noData = warpedExposure.mask.getPlaneBitMask("NO_DATA")
    warpedExposure.mask.array[~good] |= noData

    # Calibrate the variance plane
    sigNoise = warpedExposure.image.array/np.sqrt(warpedExposure.variance.array)
    rms = robustRms(sigNoise[good])
    warpedExposure.variance.array *= rms**2

    # Set the photometric calibration
    warpedExposure.setPhotoCalib(PhotoCalib(1000.0))  # 1 pixel value unit = 1 µJy = 1000 nJy

    butler.put(warpedExposure, datasetType, dataId=dataId)


def warpAsiaa(
    butler: Butler,
    exposure: ExposureF,
    timeImage: ImageF,
    band: str,
    *,
    skymapName: str = "hsc",
    datasetType: str = "deepCoadd_warp",
    sample: float = 500,
    buffer: int = 50,
    warpConfig: Warper.ConfigClass = Warper.ConfigClass(),
    workers: int = 1,
    overwrite: bool = False,
):
    """Warp an ASIAA image into our skymap system
    
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
    buffer : `int`, optional
        Number of pixels to buffer the bounding box of the warped image by when
        determining the overlap.
    warpConfig : `lsst.afw.math.Warper.ConfigClass`, optional
        Configuration for the warping.
    workers : `int`, optional
        Number of worker processes to use for warping the patches.
    overwrite : `bool`, optional
        Whether to overwrite the output if it already exists.
    """
    log = getLogger("u2k.warpAsiaa")
    registerDatasetType(butler.registry, datasetType, ("band", "tract", "patch", "skymap"), "ExposureF")
    skymap = butler.get("skyMap", collections="skymaps", skymap=skymapName)
    if exposure.metadata["BUNIT"] != "uJy":
        raise RuntimeError(f"Expected input image to be in units of uJy: got {exposure.metadata['BUNIT']}")
    wcs = exposure.getWcs()
    width, height = exposure.getDimensions()
    xNumSamples = int(width / sample)
    yNumSamples = int(height / sample)
    xSample = np.linspace(0, width - 1, xNumSamples)
    ySample = np.linspace(0, height - 1, yNumSamples)
    coordList = [wcs.pixelToSky(Point2D(x, y)) for y in ySample for x in xSample]
    patchList = defaultdict(set)
    for tract, patches in skymap.findTractPatchList(coordList):
        for pp in patches:
            patchList[tract.tract_id].add(pp.sequential_index)
    
    log.info("Found %d patches in %d tracts", sum(len(pp) for pp in patchList.values()), len(patchList))
    log.debug("Patches: %s", patchList)

    executor = ProcessPoolExecutor(workers)
    jobs = []
    for tractId in patchList:
        for patchId in patchList[tractId]:
            dataId = dict(instrument="u2k", band=band, tract=tractId, patch=patchId, skymap=skymapName)
            patchInfo = skymap[tractId][patchId]
            wcs = patchInfo.getWcs()
            bbox = Box2I()
            for vertex in patchInfo.outer_sky_polygon.getVertices():
                bbox.include(Point2I(wcs.skyToPixel(SpherePoint(vertex))))
            bbox.grow(buffer)
            bbox.clip(exposure.getBBox())
            subExposure = exposure[bbox]
            subTimeImage = timeImage[bbox]

            jobs.append(
                executor.submit(
                    warpAsiaaPatch,
                    log,
                    butler,
                    datasetType,
                    dataId,
                    subExposure,
                    subTimeImage,
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
