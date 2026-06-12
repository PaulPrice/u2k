import os

from lsst.pipe.base import Task
from lsst.pipe.tasks.calibrateImage import CalibrateImageTask
from lsst.pex.config import Config
from lsst.utils import getPackageDir

from .detector import makeDetector


class NullConfig(Config):
    """Configuration for the NullTask, which does nothing."""
    pass


class NullTask(Task):
    """A task that does nothing."""
    ConfigClass = NullConfig
    _DefaultName = "null"

    def run(self, *args, **kwargs):
        pass


class CoaddCalibrateImageConfig(CalibrateImageTask.ConfigClass):
    def setDefaults(self):
        super().setDefaults()
        self.psf_repair.retarget(NullTask)
        self.do_illumination_correction = False
        self.do_calibrate_pixels = True
        self.optional_outputs = ["stars_footprints"]
        self.run_sattle = False

        configDir = os.path.join(getPackageDir("u2k"), "config")
        self.psf_source_measurement.load(os.path.join(configDir, "apertures.py"))
        self.psf_source_measurement.load(os.path.join(configDir, "kron.py"))
        self.psf_source_measurement.load(os.path.join(configDir, "convolvedFluxes.py"))


class CoaddCalibrateImageTask(CalibrateImageTask):
    """Calibrate a coadd image"""
    ConfigClass = CoaddCalibrateImageConfig

    def run(self, exposure):
        """Calibrate a coadd image

        Parameters
        ----------
        exposure : `lsst.afw.image.ExposureF`
            The coadd image to calibrate.

        Returns
        -------
        calibrated : `lsst.afw.image.ExposureF`
            The calibrated coadd image.
        """
        exposure.setDetector(makeDetector())
        exposure.metadata["SFM_ASTROM_OFFSET_MEAN"] = 0.0
        exposure.metadata["SFM_ASTROM_OFFSET_MEDIAN"] = 0.0
        return super().run(exposures=exposure)  # Yes, plural

    def _fit_astrometry(self, *args, **kwargs):
        return [], None

    def _fit_photometry(self, exposure, stars):
        return stars, [], None, None

    def _summarize(self, *args, **kwargs):
        pass

    def _compute_per_amp_fraction(self, *args, **kwargs):
        return 0, 0.0, False
