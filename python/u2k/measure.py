from lsst.pipe.tasks.multiBand import MeasureMergedCoaddSourcesTask
from lsst.drp.tasks.forcedPhotCoadd import ForcedPhotCoaddTask


class MeasureTask(MeasureMergedCoaddSourcesTask):
    """Task to measure sources on a coadd image.
    
    This override adds required mask planes that might be missing.
    """

    def run(self, exposure, *args, **kwargs):
        maskPlanes = set()
        maskPlanes |= set(self.config.measurement.plugins["base_PixelFlags"].masksFpCenter)
        maskPlanes |= set(self.config.measurement.plugins["base_PixelFlags"].masksFpAnywhere)
        maskPlanes |= set(self.config.measurement.plugins["base_PsfFlux"].badMaskPlanes)
        for name in maskPlanes:
            exposure.mask.addMaskPlane(name)

        return super().run(exposure, *args, **kwargs)


class ForcedTask(ForcedPhotCoaddTask):
    """Task to measure sources on a coadd image.
    
    This override adds required mask planes that might be missing.
    """

    def run(self, measCat, exposure, *args, **kwargs):
        maskPlanes = set()
        maskPlanes |= set(self.config.measurement.plugins["base_PixelFlags"].masksFpCenter)
        maskPlanes |= set(self.config.measurement.plugins["base_PixelFlags"].masksFpAnywhere)
        maskPlanes |= set(self.config.measurement.plugins["base_PsfFlux"].badMaskPlanes)
        for name in maskPlanes:
            exposure.mask.addMaskPlane(name)

        return super().run(measCat, exposure, *args, **kwargs)
