from numpy.typing import ArrayLike
import astropy.table
from astro_metadata_translator.headers import merge_headers

import lsst.obs.base.utils
from lsst.utils import getPackageDir
from lsst.pipe.tasks.postprocess import (
    TransformObjectCatalogTask, TransformObjectCatalogConnections, TransformObjectCatalogConfig,
)
import lsst.pipe.tasks.postprocess

from .instrument import FILTER_PRIORITY


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
        self.outputBands = FILTER_PRIORITY


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


# The following class is not used directly, but serves as trigger to monkeypatch
# the TableVStack used by ConsolidateObjectTableTask, which is used by the pipeline.
class ConsolidateObjectTableTask(lsst.pipe.tasks.postprocess.ConsolidateObjectTableTask):
    pass


class TableVStack(lsst.obs.base.utils.TableVStack):
    def extend(self, table: astropy.table.Table, extra_values: dict[str, ArrayLike] | None = None) -> None:
        """Add a single table to the stack.

        Parameters
        ----------
        table : `astropy.table.Table`
            An astropy table instance.
        extra_values : `dict`
            Dict keyed by column name with an array-like of values to set
            for this table only.
        """
        if extra_values is None:
            extra_values = {}
        if self.result is None:
            self.result = astropy.table.Table()
            slicer = slice(None, len(table))
            for name in table.colnames:
                column = table[name]
                column_cls = type(column)
                self.result[name] = column_cls.info.new_like([column], self.capacity, name=name)
                self.result[name][: len(table)] = column
            for name, values in extra_values.items():
                self.set_extra_values(
                    table=self.result,
                    key=name,
                    values=values,
                    capacity=self.capacity,
                    slicer=slicer,
                )
            self.index = len(table)
            self.result.meta = table.meta.copy()
        else:
            next_index = self.index + len(table)
            slicer = slice(self.index, next_index)
            for name in table.colnames:
                out_col = self.result[name]
                in_col = table[name]
                if out_col.dtype != in_col.dtype:
                    # Allow string/bytes columns that differ only in their
                    # itemsize (e.g. 'U10' vs 'U12'); widen the existing
                    # column in-place if the incoming column is wider.
                    out_kind = out_col.dtype.kind
                    in_kind = in_col.dtype.kind
                    if out_kind in ("U", "S") and in_kind == out_kind:
                        if in_col.dtype.itemsize > out_col.dtype.itemsize:
                            new_dtype = in_col.dtype
                            new_col = out_col.astype(new_dtype)
                            self.result.replace_column(name, new_col)
                            out_col = self.result[name]
                        # else: out_col is already wide enough (or wider),
                        # so the incoming column will be safely cast/
                        # truncate-free on assignment.
                    else:
                        raise TypeError(
                            f"Type mismatch on column {name!r}: {out_col.dtype} != {in_col.dtype}."
                        )
                self.result[name][slicer] = table[name]
            for name, values in extra_values.items():
                self.set_extra_values(
                    table=self.result,
                    key=name,
                    values=values,
                    capacity=self.capacity,
                    slicer=slicer,
                    validate=False,
                )
            self.index = next_index
            # Butler provenance should be stripped on merge. It will be
            # added by butler on write. No attempt is made here to combine
            # provenance from multiple input tables.
            self.result.meta = merge_headers([self.result.meta, table.meta], mode="drop")
            lsst.obs.base.utils.strip_provenance_from_fits_header(self.result.meta)

lsst.obs.base.utils.TableVStack = TableVStack
lsst.pipe.tasks.postprocess.TableVStack = TableVStack
