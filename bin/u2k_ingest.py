#!/usr/bin/env python

import argparse

from lsst.daf.butler import Butler
from u2k.ingest import ingest_coadds

def main():
    parser = argparse.ArgumentParser(description="Ingest coadd files into the butler")
    parser.add_argument("butler", help="Path to the butler repository")
    parser.add_argument("path", help="Path to the coadd files")
    parser.add_argument("run", help="Run name for the ingested data")
    parser.add_argument("--datasetTypeName", default="deepCoadd", help="Dataset type name to use for the ingested data")
    parser.add_argument("--skymapName", default="hsc", help="Skymap name to use for the ingested data")
    parser.add_argument("--transfer", default="copy", help="Transfer method to use for ingesting the data")
    parser.add_argument("--id", nargs="*", action="append", help="KEY=VALUE pairs of dataId")
    parser.add_argument(
        "--skip-existing",
        default=False,
        action="store_true",
        help="Skip files that already exist in the butler?",
    )

    args = parser.parse_args()

    dataIds: list[dict[str, str|int]] | None = None
    if args.id:
        identifiers: list[list[str]] = args.id
        dataIds = [
            {key: value for key, value in [item.split("=") for item in ident]} for ident in identifiers
        ]
        for dataId in dataIds:
            for key in ("tract", "patch"):
                if key in dataId:
                    dataId[key] = int(dataId[key])

    butler = Butler(args.butler, writeable=True)
    ingest_coadds(
        butler,
        args.path,
        args.run,
        datasetTypeName=args.datasetTypeName,
        skymapName=args.skymapName,
        transfer=args.transfer,
        dataIds=dataIds,
        skipExisting=args.skip_existing,
    )


if __name__ == "__main__":
    main()
