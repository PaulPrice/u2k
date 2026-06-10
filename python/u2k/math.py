import numpy as np
from numpy.typing import ArrayLike

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
