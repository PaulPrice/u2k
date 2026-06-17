import numpy as np
from numpy.typing import ArrayLike
from lsst.afw.image import MaskedImage

__all__ = ("robustRms",)


def robustRms(array: ArrayLike, nanSafe=False) -> float:
    """Calculate a robust RMS of the array using the inter-quartile range

    Uses the standard conversion of IQR to RMS for a Gaussian.

    Parameters
    ----------
    array : `numpy.ndarray`
        Array for which to calculate RMS.
    nanSafe : `bool`, optional
        If `True`, ignore NaN values in the array.

    Returns
    -------
    rms : `float`
        Robust RMS.
    """
    if nanSafe:
        array = np.asarray(array)
        array = array[np.isfinite(array)]
    if array.size == 0:
        return np.nan
    lq, uq = np.percentile(array, (25.0, 75.0))
    return 0.741 * (uq - lq)


def renormalizeVariance(maskedImage: MaskedImage, apply: bool = True) -> float:
    """Renormalize the variance plane so that the image RMS matches expectation

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        The masked image for which to renormalize the variance.
    apply : `bool`, optional
        Apply the renormalization to the variance plane. If `False`, just
        calculate the normalization factor.

    Returns
    -------
    rms : `float`
        The square-root of the normalization factor that was applied to the
        variance plane.
    """
        
    good = np.isfinite(maskedImage.image.array)
    good &= (maskedImage.variance.array > 0)
    good &= (maskedImage.mask.array & maskedImage.mask.getPlaneBitMask("NO_DATA") == 0)

    median = np.median(maskedImage.image.array[good])
    with np.errstate(divide="ignore", invalid="ignore"):
        sigNoise = (maskedImage.image.array - median)/np.sqrt(maskedImage.variance.array)
    rms = robustRms(sigNoise[good])
    if not np.isfinite(rms) or rms <= 0:
        raise RuntimeError(f"Calculated non-positive or non-finite RMS: {rms}")
    if apply:
        maskedImage.variance.array *= rms**2
    return rms
