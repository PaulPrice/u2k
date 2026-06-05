from lsst.utils import getPackageDir
from lsst.pipe.tasks.postprocess import (
    TransformObjectCatalogTask, TransformObjectCatalogConnections, TransformObjectCatalogConfig
)


class NoEpochTransformObjectCatalogConnections(TransformObjectCatalogConnections):
    """Connections for TransformObjectCatalogTask that do not include epoch information."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        del self.inputCatalogEpoch


class NoEpochTransformObjectCatalogConfig(
    TransformObjectCatalogConfig, pipelineConnections=NoEpochTransformObjectCatalogConnections
):
    """Config for TransformObjectCatalogTask that do not include epoch information."""

    def setDefaults(self):
        super().setDefaults()
        self.goodFlags = []  # We don't have the same flags
        self.functorFile = getPackageDir("u2k") + "/schemas/Object.yaml"


class NoEpochTransformObjectCatalogTask(TransformObjectCatalogTask):
    """Task to transform an object catalog for a coadd exposure, without epoch information."""

    ConfigClass = NoEpochTransformObjectCatalogConfig
    _DefaultName = "transformObjectCatalog"
    datasets_multiband = ("ref", "Exp_multiprofit", "Sersic_multiprofit")

    def runQuantum(self, butler, inputRefs, outputRefs):
        inputs = butler.get(inputRefs)
        if self.funcs is None:
            raise ValueError("config.functorFile is None. "
                             "Must be a valid path to yaml in order to run Task as a PipelineTask.")
        result = self.run(
            handle=inputs["inputCatalog"],
            funcs=self.funcs,
            dataId=dict(outputRefs.outputCatalog.dataId.mapping),
            handle_ref=inputs["inputCatalogRef"],
            handle_Exp_multiprofit=inputs["inputCatalogExpMultiprofit"],
            handle_Sersic_multiprofit=inputs["inputCatalogSersicMultiprofit"],
        )
        butler.put(result, outputRefs)
