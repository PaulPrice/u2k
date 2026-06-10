from lsst.daf.butler import Butler, DatasetType, Registry

__all__ = ("registerDatasetType",)


def registerDatasetType(
    registry: Registry,
    name: str,
    dimensions: list[str],
    storageClass: str,
    isCalibration: bool = False,
) -> DatasetType:
    """Register a datasetType

    Parameters
    ----------
    registry : `Registry`
        Butler datastore registry.
    name : `str`
        Name of the dataset type to register.
    dimensions : list of `str`
        Relevant dimensions for the dataset type.
    storageClass : `str`
        Name of the storage class for the dataset type.
    isCalibration : `bool`, optional
        Is this dataset type a calibration?

    Returns
    -------
    datasetType : `lsst.daf.butler.DatasetType`
        Dataset type.
    """
    datasetType = DatasetType(
        name,
        dimensions,
        storageClass,
        universe=registry.dimensions,
        isCalibration=isCalibration,
    )
    registry.registerDatasetType(datasetType)
    return datasetType
