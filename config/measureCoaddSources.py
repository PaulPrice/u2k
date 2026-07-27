import os.path

config.measurement.plugins["base_PsfFlux"].badMaskPlanes = ["NO_DATA", "BAD", "SAT", "INTRP"]
config.measurement.plugins["base_PixelFlags"].masksFpCenter = [
    "CLIPPED", "SENSOR_EDGE", "INEXACT_PSF", "NO_DATA"
]
config.measurement.plugins["base_Variance"].mask = ["NO_DATA", "DETECTED", "DETECTED_NEGATIVE", "BAD", "SAT"]


config.measurement.load(os.path.join(os.path.dirname(__file__), "apertures.py"))
config.measurement.load(os.path.join(os.path.dirname(__file__), "kron.py"))
config.measurement.load(os.path.join(os.path.dirname(__file__), "convolvedFluxes.py"))
config.measurement.load(os.path.join(os.path.dirname(__file__), "hsm.py"))
config.load(os.path.join(os.path.dirname(__file__), "cmodel.py"))
