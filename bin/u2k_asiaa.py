#!/usr/bin/env python

import os
from argparse import ArgumentParser
from glob import glob

from lsst.afw.image import ImageF
from lsst.afw.math import Warper
from lsst.daf.butler import Butler

from u2k.asiaa import warpAsiaa, readAsiaa
from u2k.calibrateImage import CoaddCalibrateImageTask


def main():
    parser = ArgumentParser(description="Warp ASIAA images into our skymap system")
    parser.add_argument("butler", help="The butler path to use")
    parser.add_argument("run", help="The run to use for the warped images")
    parser.add_argument("inputDir", help="The directory containing the input ASIAA images")
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
    parser.add_argument(
        "--buffer",
        type=int,
        default=500,
        help="The number of pixels to buffer around the bounding box of the patch when warping.",
    )
    parser.add_argument("--warpConfig", help="Path to a config file for the warping.")
    parser.add_argument("--calibrateConfig", help="Path to a config file for the calibration.")
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
    calibrateConfig.do_adaptive_threshold_detection = False
    calibrateConfig.psf_detection.reEstimateBackground = True
    calibrateConfig.star_detection.reEstimateBackground = True
    if args.calibrateConfig:
        calibrateConfig.load(args.calibrateConfig)

    for pattern, band in (
        ("WI_COSMOS_?-J.fits", "WIRCAM-J"),
        ("WI_XMMLSS_?-J.fits", "WIRCAM-J"),
        ("WF_DEEP23-K-combine.fits", "WFCAM-K"),
        ("WF_ECOSMOS-??-J-combine.fits", "WFCAM-J"),
        ("WF_ECOSMOS-??-H-combine.fits", "WFCAM-H"),
        ("WF_ECOSMOS-??-K-combine.fits", "WFCAM-K"),
    ):
        mapList = sorted(glob(os.path.join(args.inputDir, "map." + pattern)))
        timeList = sorted(glob(os.path.join(args.inputDir, "time." + pattern)))
        if len(mapList) != len(timeList):
            raise RuntimeError(
                f"Found {len(mapList)} maps but {len(timeList)} time images for pattern {pattern}"
            )
        exposureList = [readAsiaa(filename) for filename in mapList]
        timeImageList = [ImageF(filename) for filename in timeList]
        warpAsiaa(
            butler,
            exposureList,
            timeImageList,
            band,
            skymapName=args.skymapName,
            datasetType=args.datasetType,
            sample=args.sample,
            buffer=args.buffer,
            calibrateConfig=calibrateConfig,
            warpConfig=warpConfig,
            workers=args.workers,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    main()
