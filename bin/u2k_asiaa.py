#!/usr/bin/env python

from argparse import ArgumentParser

from lsst.afw.image import ExposureF, ImageF
from lsst.afw.math import Warper
from lsst.daf.butler import Butler

from u2k.asiaa import warpAsiaa


def main():
    parser = ArgumentParser(description="Warp ASIAA images into our skymap system")
    parser.add_argument("butler", help="The butler path to use")
    parser.add_argument("run", help="The run to use for the warped images")
    parser.add_argument("band", help="The band of the image")
    parser.add_argument("exposure", help="The image to warp")
    parser.add_argument("timeImage", help="The time image (exposure time as a function of position). This informs the warping.")
    parser.add_argument(
        "--datasetType", default="deepCoadd_warp", help="The dataset type name for the warped image"
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
        default=50,
        help="The number of pixels to buffer around the bounding box of the patch when warping.",
    )
    parser.add_argument("--warpConfig", help="Path to a config file for the warping.")
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
    exposure = ExposureF(args.exposure)
    timeImage = ImageF(args.timeImage)

    warpConfig = Warper.ConfigClass()
    if args.warpConfig:
        warpConfig.load(args.warpConfig)

    warpAsiaa(
        butler,
        exposure,
        timeImage,
        args.band,
        skymapName=args.skymapName,
        datasetType=args.datasetType,
        sample=args.sample,
        buffer=args.buffer,
        warpConfig=warpConfig,
        workers=args.workers,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
