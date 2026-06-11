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


## Reading the data

From python:

    >>> from lsst.daf.butler import Butler
    >>> butler = Butler("/path/to/repo", collections="my_collection")
    >>> onePatch = butler.get("objectTable", skymap="hsc", tract=9812, patch=13)
    >>> fullTract = butler.get("objectTable_tract", skymap="hsc", tract=9812)

The main data products are:
* `objectTable`: a parquet file that the butler reads as a pandas `DataFrame`,
  containing all data for a patch (multiple bands merged).
* `objectTable_tract`: a similar file, but with patches within a tract merged.


## Tracts of interest

* UD-Cosmos: 9570, 9571, 9812, 9813, 9814, 10054, 10055
* UD-SXDS: 8523, 8524, 8765, 8766
* D-Cosmos: 9569, 9570, 9571, 9572, 9812, 9813, 9814, 10054, 10055, 10056
* D-XMM-LSS: 8282, 8283, 8284, 8523, 8524, 8525, 8765, 8766, 8767
* D-ELAIS-N1: 16984, 16985, 17129, 17130, 17131, 17270, 17271, 17272, 17406, 17407
* D-DEEP2-3: 9219, 9220, 9221, 9462, 9463, 9464, 9465, 9706, 9707, 9708


## New data

### ASIAA

ASIAA have provided coadded images from CFHT-WIRCAM and UKIRT-WFCAM,
along with exposure time maps.

    u2k_asiaa.py $REPO asiaa WFCAM-K map.WF_DEEP23-K-combine.fits time.WF_DEEP23-K-combine.fits -j 10 > $REPO/ingest/WF_DEEP23-K.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-H map.WF_ECOSMOS-NE-H-combine.fits time.WF_ECOSMOS-NE-H-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-NE-H.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-K map.WF_ECOSMOS-NE-K-combine.fits time.WF_ECOSMOS-NE-K-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-NE-K.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-H map.WF_ECOSMOS-NW-H-combine.fits time.WF_ECOSMOS-NW-H-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-NW-H.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-J map.WF_ECOSMOS-NW-J-combine.fits time.WF_ECOSMOS-NW-J-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-NW-J.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-K map.WF_ECOSMOS-NW-K-combine.fits time.WF_ECOSMOS-NW-K-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-NW-K.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-H map.WF_ECOSMOS-SE-H-combine.fits time.WF_ECOSMOS-SE-H-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-SE-H.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-J map.WF_ECOSMOS-SE-J-combine.fits time.WF_ECOSMOS-SE-J-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-SE-J.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-K map.WF_ECOSMOS-SE-K-combine.fits time.WF_ECOSMOS-SE-K-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-SE-K.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-H map.WF_ECOSMOS-SW-H-combine.fits time.WF_ECOSMOS-SW-H-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-SW-H.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-J map.WF_ECOSMOS-SW-J-combine.fits time.WF_ECOSMOS-SW-J-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-SW-J.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-K map.WF_ECOSMOS-SW-K-combine.fits time.WF_ECOSMOS-SW-K-combine.fits -j 10 > $REPO/ingest/WF_ECOSMOS-SW-K.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-WIRCAM-J map.WF+WI_COSMOS-J.fits time.WF+WI_COSMOS-J.fits -j 10 > $REPO/ingest/WF+WI_COSMOS-J.log 2>&1
    u2k_asiaa.py $REPO asiaa WFCAM-WIRCAM-J map.WF+WI_DEEP23-J.fits time.WF+WI_DEEP23-J.fits -j 10 > $REPO/ingest/WF+WI_DEEP23-J.log 2>&1
    u2k_asiaa.py $REPO asiaa WIRCAM-J map.WI_COSMOS_a-J.fits time.WI_COSMOS_a-J.fits -j 10 > $REPO/ingest/WI_COSMOS_a-J.log 2>&1
    u2k_asiaa.py $REPO asiaa WIRCAM-J map.WI_COSMOS_d-J.fits time.WI_COSMOS_d-J.fits -j 10 > $REPO/ingest/WI_COSMOS_d-J.log 2>&1
    u2k_asiaa.py $REPO asiaa WIRCAM-J map.WI_XMMLSS_N-J.fits time.WI_XMMLSS_N-J.fits -j 10 > $REPO/ingest/WI_XMMLSS_N-J.log 2>&1
    u2k_asiaa.py $REPO asiaa WIRCAM-J map.WI_XMMLSS_S-J.fits time.WI_XMMLSS_S-J.fits -j 10 > $REPO/ingest/WI_XMMLSS_S-J.log 2>&1

### VISTA VEILS

The ESO data is delivered with munged filenames, and need to get renamed:

    cat readme_56f96d63-b148-48c3-8591-b5df83c13508.txt | awk '$2 ~ /^ADP/ {print "test -f",$2,"&& mv",$2,$3}' | sh

Then ingest:

    u2k_veils.py $REPO veils . -j 10

(The code iterates over the bands, using hardwired globs.)
