from collections import defaultdict

from lsst.pipe.tasks.postprocess import (
    WriteObjectTableConnections, WriteObjectTableTask, WriteObjectTableConfig
)
from lsst.pipe.base import Struct


class NoMultiprofitWriteObjectTableConnections(WriteObjectTableConnections):
    """Connections for NoMultiprofitWriteObjectTableTask"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        del self.inputCatalogPsfsMultiprofit


class NoMultiprofitWriteObjectTableConfig(
    WriteObjectTableConfig, pipelineConnections=NoMultiprofitWriteObjectTableConnections
):
    """Config for NoMultiprofitWriteObjectTableTask"""
    pass


class NoMultiprofitWriteObjectTableTask(WriteObjectTableTask):
    """Task to write an object table for a coadd exposure, without multiprofit catalogs."""

    ConfigClass = NoMultiprofitWriteObjectTableConfig

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        catalogs = defaultdict(dict)
        for dataset, connection in (
            ("meas", "inputCatalogMeas"),
            ("forced_src", "inputCatalogForcedSrc"),
        ):
            for ref, cat in zip(getattr(inputRefs, connection), inputs[connection]):
                catalogs[ref.dataId["band"]][dataset] = cat

        dataId = butlerQC.quantum.dataId
        df = self.run(catalogs=catalogs, tract=dataId["tract"], patch=dataId["patch"])
        outputs = Struct(outputCatalog=df)
        butlerQC.put(outputs, outputRefs)
