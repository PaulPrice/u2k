# Enable HSM shapes (unsetup meas_extensions_shapeHSM to disable)
# 'config' is a SourceMeasurementConfig.
import os.path
from lsst.utils import getPackageDir

try:
    config.load(os.path.join(getPackageDir("meas_extensions_shapeHSM"), "config", "enable.py"))
    config.plugins["ext_shapeHSM_HsmShapeRegauss"].deblendNChild = "deblend_nChild"
    # Enable debiased moments
    config.plugins.names |= ["ext_shapeHSM_HsmPsfMomentsDebiased"]
except LookupError as e:
    print("Cannot enable shapeHSM (%s): disabling HSM shape measurements" % (e,))

for name in (
    "ext_shapeHSM_HigherOrderMomentsSource",
    "ext_shapeHSM_HsmSourceMoments",
    "ext_shapeHSM_HsmSourceMomentsRound",
    "ext_shapeHSM_HsmPsfMomentsDebiased",
    "ext_shapeHSM_HsmShapeBj",
    "ext_shapeHSM_HsmShapeLinear",
    "ext_shapeHSM_HsmShapeKsb",
    "ext_shapeHSM_HsmShapeRegauss",
):
    config.plugins[name].badMaskPlanes = ["NO_DATA", "BAD", "SAT", "INTRP"]
