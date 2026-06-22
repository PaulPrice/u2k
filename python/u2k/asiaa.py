from logging import Logger

import numpy as np

from lsst.afw.image import ExposureF, ImageF, PhotoCalib, FilterLabel
from lsst.afw.math import Warper
from lsst.daf.butler import Butler
from lsst.geom import Point2I, Box2I, SpherePoint
from lsst.skymap import PatchInfo

from .butler import registerDatasetType
from .calibrateImage import CoaddCalibrateImageTask
from .logging import getLogger
from .math import renormalizeVariance
from .processPool import ProcessPool
from .warp import getWarpMapping, warpExposures


def readAsiaa(filename: str) -> ExposureF:
    """Read an ASIAA image.

    There's something wrong that happens when the LSST afw Exposure code reads
    the ASIAA coadd images directly.
    
    Parameters
    ----------
    filename : `str`
        The path to the image file.
    """
    import astropy.io.fits
    from lsst.afw.fits import readMetadata
    from lsst.afw.image import makeMaskedImage, makeExposure
    from lsst.afw.geom import SkyWcs
    with astropy.io.fits.open(filename) as fits:
        image = fits[0].data.astype(np.float32)
    metadata = readMetadata(filename)
    wcs = SkyWcs(metadata)
    image = ImageF(image)
    exposure = makeExposure(makeMaskedImage(image), wcs)
    exposure.setMetadata(metadata)
    return exposure


def calibrateAsiaa(
    log: Logger,
    exposure: ExposureF,
    timeImage: ImageF,
    calibrateConfig: CoaddCalibrateImageTask.ConfigClass,
):
    """Calibrate an ASIAA image
    
    Parameters
    ----------
    log : `logging.Logger`
        Logger to use for logging messages.
    exposure : `lsst.afw.image.ExposureF`
        The image to calibrate.
    timeImage : `lsst.afw.image.ImageF`
        The time image (exposure time as a function of position). This informs the
        variance plane.
    calibrateConfig : `u2k.CoaddCalibrateImageTask.ConfigClass`
        Configuration for the calibration.
    """
    exposure.mask.array[:] = 0
    # The variance goes as the background, which goes as the exposure time
    # Somehow, some of the times can be negative...
    exposure.variance.array[:] = np.clip(timeImage.array, 0.0, None)

    good = np.isfinite(exposure.image.array)
    exposure.image.array[~good] = 0.0
    good &= np.isfinite(exposure.variance.array) & (exposure.variance.array > 0)

    noData = exposure.mask.getPlaneBitMask("NO_DATA")
    exposure.mask.array[~good] |= noData
    assert(np.all(exposure.variance.array >= 0))

    rms = renormalizeVariance(exposure.maskedImage)
    log.info("Renormalized variance plane with RMS factor %.3f", rms)

    task = CoaddCalibrateImageTask(config=calibrateConfig, log=log)
    task.run(exposure=exposure)


def warpAsiaaPatch(
    log: Logger,
    butler: Butler,
    datasetType: str,
    dataId: dict[str, str|int],
    exposureList: list[ExposureF],
    timeImageList: list[ImageF],
    patchInfo: PatchInfo,
    warpConfig: Warper.ConfigClass,
    calibrateConfig: CoaddCalibrateImageTask.ConfigClass,
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
    exposureList : `list[lsst.afw.image.ExposureF]`
        The images to warp.
    timeImageList : `list[lsst.afw.image.ImageF]`
        The time images (exposure times as a function of position). This informs
        the variance plane.
    patchInfo : `lsst.skymap.PatchInfo`
        The patch to warp into.
    warpConfig : `lsst.afw.math.Warper.ConfigClass`
        Configuration for the warping.
    calibrateConfig : `u2k.CoaddCalibrateImageTask.ConfigClass`
        Configuration for the calibration.
    overwrite : `bool`, optional
        Whether to overwrite the output if it already exists.
    """
    assert len(exposureList) == len(timeImageList), "exposureList and timeImageList must have the same length"

    dataRef = butler.find_dataset(datasetType, dataId)
    if dataRef:
        if overwrite:
            log.info("Overwriting existing output for %s", dataId)
            butler.pruneDatasets([dataRef], disassociate=True, unstore=True, purge=True)
        else:
            log.info("Output for %s already exists and overwrite=False, skipping", dataId)
            return

    for ii, (exposure, timeImage) in enumerate(zip(exposureList, timeImageList)):
        log.info("Calibrating image %d for patch %s", ii, dataId)
        try:
            calibrateAsiaa(log, exposure, timeImage, calibrateConfig)
        except Exception as exc:
            log.error("Failed to calibrate image %d for patch %s: %s", ii, dataId, exc)
            return

    log.info("Warping %d exposures into patch %s", len(exposureList), dataId)
    warped = warpExposures(
        log=log,
        exposureList=exposureList,
        patchInfo=patchInfo,
        warpConfig=warpConfig,
        datasetType=datasetType,
        dataId=dataId,
    )
    if warped is not None:
        butler.put(warped, datasetType, dataId=dataId)


def warpAsiaa(
    butler: Butler,
    exposureList: list[ExposureF],
    timeImageList: list[ImageF],
    band: str,
    *,
    skymapName: str = "hsc",
    datasetType: str = "deepCoadd",
    sample: float = 500,
    buffer: int = 500,
    warpConfig: Warper.ConfigClass = Warper.ConfigClass(),
    calibrateConfig: CoaddCalibrateImageTask.ConfigClass = CoaddCalibrateImageTask.ConfigClass(),
    workers: int = 1,
    overwrite: bool = False,
):
    """Warp an ASIAA image into our skymap system
    
    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        The butler to use.
    exposureList : `list[lsst.afw.image.ExposureF]`
        The images to warp.
    timeImageList : `list[lsst.afw.image.ImageF]`
        The time images (exposure time as a function of position). This informs
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
    calibrateConfig : `u2k.CoaddCalibrateImageTask.ConfigClass`, optional
        Configuration for the calibration.
    workers : `int`, optional
        Number of worker processes to use for warping the patches.
    overwrite : `bool`, optional
        Whether to overwrite the output if it already exists.
    """
    log = getLogger("u2k.warpAsiaa")
    registerDatasetType(butler.registry, datasetType, ("band", "tract", "patch", "skymap"), "ExposureF")
    skymap = butler.get("skyMap", collections="skymaps", skymap=skymapName)

    for exposure, timeImage in zip(exposureList, timeImageList):
        if exposure.metadata["BUNIT"] != "uJy":
            raise RuntimeError(
                f"Expected input image to be in units of uJy: got {exposure.metadata['BUNIT']}"
            )
        exposure.setPhotoCalib(PhotoCalib(1000.0))  # 1 ADU = 1 µJy = 1000 nJy
        exposure.setFilter(FilterLabel.fromBandPhysical(band, band))

        exposure.variance.array[:] = np.clip(timeImage.array, 0.0, None)
        exposure.mask[:] = 0

    patchToIndices = getWarpMapping(
        [exposure.getWcs() for exposure in exposureList],
        [exposure.getDimensions() for exposure in exposureList],
        skymap,
        sample=sample,
        log=log,
    )    

    pool = ProcessPool(workers)
    with pool:
        for tractId, patchId in patchToIndices:
            dataId = dict(instrument="u2k", band=band, tract=tractId, patch=patchId, skymap=skymapName)
            patchInfo = skymap[tractId][patchId]

            subExposureList = []
            subTimeImageList = []
            for ii in patchToIndices[(tractId, patchId)]:
                exposure = exposureList[ii]
                timeImage = timeImageList[ii]
                bbox = Box2I()  # Bounding box of patch on the exposure, in pixel coordinates
                for vertex in patchInfo.outer_sky_polygon.getVertices():
                    bbox.include(Point2I(exposure.getWcs().skyToPixel(SpherePoint(vertex))))
                bbox.clip(exposure.getBBox())
                if bbox.isEmpty():
                    continue
                bbox.grow(buffer)  # Extra pixels to get a good PSF model
                bbox.clip(exposure.getBBox())
                if bbox.isEmpty():
                    continue

                subExposureList.append(exposure[bbox].clone())
                subTimeImageList.append(timeImage[bbox])
            if not subExposureList:
                continue

            pool.submit(
                warpAsiaaPatch,
                log,
                butler,
                datasetType,
                dataId,
                subExposureList,
                subTimeImageList,
                patchInfo,
                warpConfig,
                calibrateConfig,
                overwrite,
            )
