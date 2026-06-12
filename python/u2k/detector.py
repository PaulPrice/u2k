from lsst.afw.cameraGeom import Camera


def makeCamera():
    """Make a camera with a single detector"""
    builder = Camera.Builder("camera")
    builder.add("detector", 0)
    camera = builder.finish()
    return camera


def makeDetector():
    """Return a detector"""
    return makeCamera()[0]
