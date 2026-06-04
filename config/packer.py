# Band listing for ID packing
from u2k.instrument import FILTER_DEFINITIONS
config.idGenerator.packer.bands = {band.band: ii for ii, band in enumerate(FILTER_DEFINITIONS)}
