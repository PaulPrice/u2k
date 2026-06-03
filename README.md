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

butler create $REPO
butler register-instrument $REPO u2k.instrument.GenericCoadd
butler register-skymap $REPO -C $U2K_DIR/config/skymap_rings.py -c name=hsc
u2k_ingest.py $REPO /path/to/data coadds

