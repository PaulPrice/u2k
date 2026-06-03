import os
import re

from lsst.daf.butler import Butler, DatasetRef, FileDataset
from .logging import getLogger


def ingest_coadds(
    butler: Butler,
    path: str,
    run: str,
    datasetTypeName: str = "deepCoadd",
    skymapName: str = "hsc",
    transfer: str = "copy",
    dataIds: list[dict[str, str|int]] | None = None,
    skipExisting: bool = False,
) -> None:
    """Ingest coadd files into the butler

    The coadd files come from previous Gen2 processing, and have filenames of
    the form ``calexp-{band}-{tract}-{xPatch},{yPatch}.fits``.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        The butler into which to ingest the data.
    path : `str`
        Path to the coadd files.
    run : `str`
        Run name for the ingested data.
    datasetTypeName : `str`, optional
        Dataset type name to use for the ingested data.
    skymapName : `str`, optional
        Skymap name to use for the ingested data.
    transfer : `str`, optional
        Transfer method to use for ingesting the data.
    dataIds : `list` of `dict`, optional
        Data IDs to restrict the data being ingested.
    skipExisting : `bool`, optional
        Skip files that already exist in the butler?
    """
    log = getLogger("u2k.ingest")
    skymap = butler.get("skyMap", collections="skymaps", skymap=skymapName)

    datasetType = butler.get_dataset_type(datasetTypeName)

    datasets = []
    for dirpath, _, filenames in os.walk(path):
        for filename in filenames:
            if not filename.startswith("calexp-") or not filename.endswith(".fits"):
                continue
            log.debug("Processing file %s", os.path.join(dirpath, filename))
            match = re.match(r"calexp-(.*)-(\d+)-(\d+),(\d+).fits", filename)
            if not match:
                continue
            band, tract, xPatch, yPatch = match.groups()
            tract = int(tract)
            xPatch = int(xPatch)
            yPatch = int(yPatch)
            patchNum = skymap[tract][(xPatch, yPatch)].sequential_index

            if band.startswith("NB"):
                band = "HSC-" + band

            ingestId = dict(instrument="u2k", band=band, tract=tract, patch=patchNum, skymap=skymapName)
            if dataIds is not None:
                found = False
                for ident in dataIds:
                    if all(ingestId[key] == ident[key] for key in ident):
                        found = True
                        break
                if not found:
                    log.info(f"Skipping {filename} because it does not match any specified dataId")
                    continue

            ref = DatasetRef(datasetType, ingestId, run)
            datasets.append(FileDataset(os.path.join(dirpath, filename), refs=[ref]))

    log.info("Ingesting %d datasets", len(datasets))
    butler.ingest(*datasets, transfer=transfer, skip_existing=skipExisting)
    log.info("Ingest complete")
