from collections.abc import Iterable
from typing import cast

from lsst.daf.butler import DeferredDatasetHandle
import lsst.analysis.tools.tasks
from lsst.analysis.tools.interfaces import KeyedData


class ObjectTableTractAnalysisTask(lsst.analysis.tools.tasks.ObjectTableTractAnalysisTask):
    def loadData(self, handle: DeferredDatasetHandle, names: Iterable[str] | None = None) -> KeyedData:
        """Load the minimal set of keyed data from the input dataset.

        Parameters
        ----------
        handle : `DeferredDatasetHandle`
            Handle to load the dataset with only the specified columns.
        names : `Iterable` of `str`
            The names of keys to extract from the dataset.
            If `names` is `None` then the `collectInputNames` method
            is called to generate the names.
            For most purposes these are the names of columns to load from
            a catalog or data frame.

        Returns
        -------
        result: `KeyedData`
            The dataset with only the specified keys loaded.
        """
        if names is None:
            names = self.collectInputNames()
        have = set(handle.get(component="columns"))
        missing = set(names) - have
        if missing:
            self.log.warning("Dropping unavailable columns: %s", ", ".join(missing))
            names = have.intersection(names)
        return cast(KeyedData, handle.get(parameters={"columns": names}))
