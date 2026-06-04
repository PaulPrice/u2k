import os

from lsst.utils import getPackageDir
from lsst.daf.butler import DatasetType, Registry
from lsst.obs.base import Instrument, FilterDefinition, FilterDefinitionCollection
from lsst.obs.base.formatters.fitsExposure import FitsExposureFormatter
from lsst.utils.introspection import get_full_type_name

__all__ = ("GenericCoadd", "FILTER_DEFINITIONS", "FILTER_PRIORITY")


FILTER_DEFINITIONS = FilterDefinitionCollection(
    FilterDefinition("HSC-G", "HSC-G", "Hyper Suprime-Cam g-band"),
    FilterDefinition("HSC-R", "HSC-R", "Hyper Suprime-Cam r-band"),
    FilterDefinition("HSC-I", "HSC-I", "Hyper Suprime-Cam I-band"),
    FilterDefinition("HSC-Z", "HSC-Z", "Hyper Suprime-Cam z-band"),
    FilterDefinition("HSC-Y", "HSC-Y", "Hyper Suprime-Cam Y-band"),
    FilterDefinition("HSC-NB0387", "HSC-NB0387", "Hyper Suprime-Cam NB0387 narrowband"),
    FilterDefinition("HSC-NB0527", "HSC-NB0527", "Hyper Suprime-Cam NB0527 narrowband"),
    FilterDefinition("HSC-NB0718", "HSC-NB0718", "Hyper Suprime-Cam NB0718 narrowband"),
    FilterDefinition("HSC-NB0816", "HSC-NB0816", "Hyper Suprime-Cam NB0816 narrowband"),
    FilterDefinition("HSC-NB0921", "HSC-NB0921", "Hyper Suprime-Cam NB0921 narrowband"),
    FilterDefinition("HSC-NB0973", "HSC-NB0973", "Hyper Suprime-Cam NB0973 narrowband"),
    FilterDefinition("HSC-NB1010", "HSC-NB1010", "Hyper Suprime-Cam NB1010 narrowband"),
    FilterDefinition("VIRCAM-H", "VIRCAM-H", "VISTA Infrared Camera H-band"),
    FilterDefinition("VIRCAM-J", "VIRCAM-J", "VISTA Infrared Camera J-band"),
    FilterDefinition("VIRCAM-Ks", "VIRCAM-Ks", "VISTA Infrared Camera Ks-band"),
    FilterDefinition("VIRCAM-NB118", "VIRCAM-NB118", "VISTA Infrared Camera NB118 narrowband"),
    FilterDefinition("VIRCAM-Y", "VIRCAM-Y", "VISTA Infrared Camera Y-band"),
    FilterDefinition("WFCAM-H", "WFCAM-H", "UKIRT Wide Field Camera H-band"),
    FilterDefinition("WFCAM-J", "WFCAM-J", "UKIRT Wide Field Camera J-band"),
    FilterDefinition("WFCAM-K", "WFCAM-K", "UKIRT Wide Field Camera K-band"),
    FilterDefinition("MegaCam-u", "MegaCam-u", "CFHT MegaCam u-band"),
    FilterDefinition("MegaCam-uS", "MegaCam-uS", "CFHT MegaCam u*-band"),
)


FILTER_PRIORITY = [
    "HSC-I", "HSC-R", "HSC-Z", "HSC-Y", "HSC-G",
    "VIRCAM-Y", "VIRCAM-J", "WFCAM-J", "VIRCAM-H", "WFCAM-H", "VIRCAM-Ks", "WFCAM-K",
    "MegaCam-uS", "MegaCam-u",
    "HSC-NB0816", "HSC-NB0718", "HSC-NB0527", "HSC-NB0921", "HSC-NB0973", "HSC-NB0387", "HSC-NB1010",
    "VIRCAM-NB118",
]


class GenericCoadd(Instrument):
    """Gen3 Butler specialization class for Subaru's Prime Focus Spectrograph."""

    policyName = "u2k"
    obsDataPackage = "u2k"
    filterDefinitions = FILTER_DEFINITIONS

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        packageDir = getPackageDir("u2k")
        self.configPaths = [os.path.join(packageDir, "config")]

    @classmethod
    def getName(cls) -> str:
        # Docstring inherited from Instrument.getName
        return "U2K"

    def writeCameraGeom(self, *args, **kwargs):
        return

    def writeStandardTextCuratedCalibrations(self, *args, **kwargs):
        return

    def getRawFormatter(self, dataId):
        return FitsExposureFormatter

    def register(self, registry: Registry, update: bool = False) -> None:
        with registry.transaction():
            registry.syncDimensionData(
                "instrument",
                {
                    "name": "U2K",
                    "detector_max": 0,
                    "visit_max": 0,
                    "class_name": get_full_type_name(self),
                },
                update=update,
            )
            self._registerFilters(registry, update=update)
        self.registerDatasetType(registry, "deepCoadd", ["band", "tract", "patch", "skymap"], "ExposureF")

    def registerDatasetType(
        self,
        registry: Registry,
        name: str,
        dimensions: list[str],
        storageClass: str,
        isCalibration: bool = False,
    ):
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

    def getCamera(self):
        raise NotImplementedError("We don't have a camera definition")
