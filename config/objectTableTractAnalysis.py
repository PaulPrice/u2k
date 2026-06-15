from lsst.analysis.tools.actions.vector import ConvertFluxToMag
from lsst.analysis.tools.atools import WPerpCModel, WPerpPSF, XPerpCModel, XPerpPSF, YPerpCModel, YPerpPSF

for name, cls, bands, fluxType in (
    ("wPerpPSF", WPerpPSF, ("HSC-G", "HSC-R", "HSC-I"), "psfFlux"),
    ("wPerpCModel", WPerpCModel, ("HSC-G", "HSC-R", "HSC-I"), "cModelFlux"),
    ("xPerpPSF", XPerpPSF, ("HSC-G", "HSC-R", "HSC-I"), "psfFlux"),
    ("xPerpCModel", XPerpCModel, ("HSC-G", "HSC-R", "HSC-I"), "cModelFlux"),
    ("yPerpPSF", YPerpPSF, ("HSC-R", "HSC-I", "HSC-Z"), "psfFlux"),
    ("yPerpCModel", YPerpCModel, ("HSC-R", "HSC-I", "HSC-Z"), "cModelFlux"),
):
    setattr(config.atools, name, cls)
    conf = getattr(config.atools, name)
    conf.prep.selectors.flagSelector.bands = bands
    conf.prep.selectors.snSelector.bands = [bands[1]]
    conf.prep.selectors.magSelector.bands = [bands[1]]
    conf.prep.selectors.starSelector1.vectorKey = bands[0] + "_extendedness"
    conf.prep.selectors.starSelector2.vectorKey = bands[1] + "_extendedness"
    conf.prep.selectors.starSelector3.vectorKey = bands[2] + "_extendedness"
    conf.process.buildActions.x.magDiff.col1 = bands[0] + "_" + fluxType
    conf.process.buildActions.x.magDiff.col2 = bands[1] + "_" + fluxType
    conf.process.buildActions.y.magDiff.col1 = bands[1] + "_" + fluxType
    conf.process.buildActions.y.magDiff.col2 = bands[2] + "_" + fluxType
    conf.process.buildActions.mag = ConvertFluxToMag(vectorKey=bands[1] + "_" + fluxType)
