from lsst.pipe.tasks.multiBand import MeasureMergedCoaddSourcesTask
from lsst.drp.tasks.forcedPhotCoadd import ForcedPhotCoaddTask
from lsst.meas.algorithms import CoaddBoundedField


def checkExposure(exposure, wcs, plugins):
    """Check that the exposure is valid for measurement

    Fix it if not.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to check.
    wcs : `lsst.afw.geom.Wcs`
        World coordinate system to use.
    plugins : `lsst.meas.base.Plugins`
        Measurement plugins being used.
    """
    # Ensure we have the required mask planes
    maskPlanes = set()
    maskPlanes |= set(plugins["base_PixelFlags"].masksFpCenter)
    maskPlanes |= set(plugins["base_PixelFlags"].masksFpAnywhere)
    maskPlanes |= set(plugins["base_PsfFlux"].badMaskPlanes)
    for name in maskPlanes:
        exposure.mask.addMaskPlane(name)

    # Ensure we have a WCS
    if exposure.getWcs() is None:
        exposure.setWcs(wcs)

    # Ensure we have a coadd apcorr map with the correct WCS
    acm = exposure.getInfo().getApCorrMap()
    for key in acm:
        field = acm[key]
        if not isinstance(field, CoaddBoundedField) or field.getCoaddWcs() is not None:
            continue
        acm[key] = CoaddBoundedField(field.getBBox(), wcs, field.getElements(), field.getDefault())


class MeasureTask(MeasureMergedCoaddSourcesTask):
    """Task to measure sources on a coadd image.
    
    This override adds required mask planes that might be missing.
    """

    def run(self, exposure, sources, skyInfo, *args, **kwargs):
        checkExposure(exposure, skyInfo.wcs, self.config.measurement.plugins)
        return super().run(exposure, sources, skyInfo, *args, **kwargs)


class ForcedTask(ForcedPhotCoaddTask):
    """Task to measure sources on a coadd image.
    
    This override adds required mask planes that might be missing.
    """

    def run(self, measCat, exposure, refCat, refWcs, *args, **kwargs):
        checkExposure(exposure, refWcs, self.config.measurement.plugins)
        return super().run(measCat, exposure, refCat, refWcs, *args, **kwargs)
