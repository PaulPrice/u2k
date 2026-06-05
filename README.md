# Processing U2K version 3

U2K is a deep imaging dataset spanning the optical and NIR,
from u-band to K-band (hence the name),
centered around the Hyper Suprime-Cam SSP deep and ultra-deep fields:
Cosmos (D/UD), SXDS (UD), XMM-LSS (D), ELAIS-N1 (D) and DEEP2-3 (D).

This repository contains scripts and notes in support of the version 3 reprocessing of these data.
Version 1 used an early version of hscPipe, which had problems with background and source detection.
Version 2 used the last version of hscPipe (v8), but did not include some data that is now available.
Our reprocessing will utilize the latest version of the Rubin/LSST multiband pipeline,
and include new data.


## Data files

Data that were previously processed as part of version 2 are laid out in a directory structure
specified by the old LSST Gen2 butler; we need to ingest these into a modern Gen3 butler.
New data is provided as coadd FITS files constructed with an external pipeline;
these need to be reformatted (warped to the skymap, PSF and variance measured) and ingested.


## Butler setup

### Small-scale butler with SQLite

This butler is not suitable for large-scale, parallel processing runs,
but may be helpful for testing on a single tract.

    butler create $REPO
    butler register-instrument $REPO u2k.instrument.GenericCoadd
    butler register-skymap $REPO -C $U2K_DIR/config/skymap_rings.py -c name=hsc
    u2k_ingest.py $REPO /path/to/data coadds --id tract=9812  # only ingesting a single tract!

### Production butler with PostgreSQL

Set up `~/.lsst/db-auth.yaml` and `~/.pgpass`

    butler create $REPO --seed-config $U2K_DIR/butler/butler.yaml --override
    butler register-instrument $REPO u2k.instrument.GenericCoadd
    butler register-skymap $REPO -C $U2K_DIR/config/skymap_rings.py -c name=hsc
    u2k_ingest.py $REPO /path/to/data coadds


## Running the pipeline

    pipetask run -b $REPO -i coadds,skymaps -o RERUN_NAME -p '$U2K_DIR/pipelines/multiband.yaml' -d "tract = 12345 AND patch = 67" --register-dataset-types

You can drop the `--register-dataset-types` after running the pipeline once.
