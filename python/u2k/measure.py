from lsst.pipe.tasks.multiBand import MeasureMergedCoaddSourcesTask
from lsst.drp.tasks.forcedPhotCoadd import ForcedPhotCoaddTask


class MeasureTask(MeasureMergedCoaddSourcesTask):
    """Task to measure sources on a coadd image.
    
    This override adds required mask planes that might be missing.
    """

    def run(self, exposure, sources, skyInfo, *args, **kwargs):
        maskPlanes = set()
        maskPlanes |= set(self.config.measurement.plugins["base_PixelFlags"].masksFpCenter)
        maskPlanes |= set(self.config.measurement.plugins["base_PixelFlags"].masksFpAnywhere)
        maskPlanes |= set(self.config.measurement.plugins["base_PsfFlux"].badMaskPlanes)
        for name in maskPlanes:
            exposure.mask.addMaskPlane(name)

        exposure.setWcs(skyInfo.wcs)
        return super().run(exposure, sources, skyInfo, *args, **kwargs)


class ForcedTask(ForcedPhotCoaddTask):
    """Task to measure sources on a coadd image.
    
    This override adds required mask planes that might be missing.
    """

    def run(self, measCat, exposure, refCat, refWcs, *args, **kwargs):
        maskPlanes = set()
        maskPlanes |= set(self.config.measurement.plugins["base_PixelFlags"].masksFpCenter)
        maskPlanes |= set(self.config.measurement.plugins["base_PixelFlags"].masksFpAnywhere)
        maskPlanes |= set(self.config.measurement.plugins["base_PsfFlux"].badMaskPlanes)
        for name in maskPlanes:
            exposure.mask.addMaskPlane(name)

        exposure.setWcs(refWcs)
        return super().run(measCat, exposure, refCat, refWcs, *args, **kwargs)
