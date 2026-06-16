#!/usr/bin/env python

import os
from argparse import ArgumentParser
from glob import glob

from lsst.afw.math import Warper
from lsst.daf.butler import Butler

from u2k.veils import warpVeils, readVeils, readVeilsConf
from u2k.calibrateImage import CoaddCalibrateImageTask


def main():
    parser = ArgumentParser(description="Warp ASIAA images into our skymap system")
    parser.add_argument("butler", help="The butler path to use")
    parser.add_argument("run", help="The run to use for the warped images")
    parser.add_argument("source", help="Source directory for the images")
    parser.add_argument(
        "--datasetType", default="deepCoadd", help="The dataset type name for the warped image"
    )
    parser.add_argument("--skymapName", default="hsc", help="The name of the skymap to warp into")
    parser.add_argument(
        "--sample",
        type=float,
        default=500,
        help="The spacing at which to sample the input image when determining the overlap.",
    )
    parser.add_argument("--warpConfig", help="Path to a config file for the warping.")
    parser.add_argument("--calibrateConfig", help="Path to a config file for the calibration.")
    parser.add_argument(
        "--oldCollection",
        default="coadds",
        help="Collection to look in for existing coadds to combine with.",
    )
    parser.add_argument(
        "--workers",
        "-j",
        type=int,
        default=1,
        help="Number of worker processes to use for warping the patches.",
    )
    parser.add_argument(
        "--overwrite",
        default=False,
        action="store_true",
        help="Overwrite existing warped images?",
    )
    args = parser.parse_args()

    butler = Butler(args.butler, run=args.run, collections=args.run, writeable=True)

    warpConfig = Warper.ConfigClass()
    if args.warpConfig:
        warpConfig.load(args.warpConfig)
    calibrateConfig = CoaddCalibrateImageTask.ConfigClass()
    if args.calibrateConfig:
        calibrateConfig.load(args.calibrateConfig)

    for band, expGlob, confGlob in (
        ("VIRCAM-J", "XMM_Pointing_?_3?x_J_p?.fits.fz", "XMM_Pointing_?_3?x_J_p?_conf.fits.fz"),
        ("VIRCAM-Ks", "XMM_Pointing_?_3?x_Ks_p?.fits.fz", "XMM_Pointing_?_3?x_Ks_p?_conf.fits.fz"),
    ):
        exposureNames = sorted(glob(os.path.join(args.source, expGlob)))
        confNames = sorted(glob(os.path.join(args.source, confGlob)))
        if len(exposureNames) != len(confNames):
            raise ValueError(
                f"Number of exposures ({len(exposureNames)}) and confs ({len(confNames)}) "
                f"do not match for band {band}"
            )
        exposureList = []
        confList = []
        for expName, confName in zip(exposureNames, confNames):
            exposureList += readVeils(expName)
            confList += readVeilsConf(confName)

        warpVeils(
            butler,
            exposureList,
            confList,
            band,
            skymapName=args.skymapName,
            datasetType=args.datasetType,
            sample=args.sample,
            warpConfig=warpConfig,
            calibrateConfig=calibrateConfig,
            workers=args.workers,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    main()
