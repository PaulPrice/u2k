from lsst.afw.image import ExposureF, PhotoCalib

__all__ = ("fixCoadd",)


def fixCoadd(exposure: ExposureF) -> None:
    """Fix a coadd exposure to be in nJy and have no DETECTED pixels

    We're operating on ``deepCoadd_calexp`` exposures from the U2K version 2
    processing; these have already had detection run on them. The old
    ``DETECTED`` mask plane is getting in the way of our new detection, so
    we clear it first.

    We also need to rescale the image so that it's in units of nJy (not with
    a magnitude zero point of 27 as previous LSST versions used), as this
    is what the TransformObjectCatalogTask expects.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The coadd exposure to fix
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
