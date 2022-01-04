"""
This package provides helper routines for the SUVI (Solar UltraViolet Imager) 
instrument on the GOES-R series of satellites.

.. note::
    SUVI data (and data from other GOES instrumentation) can be found 
    `here <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/>`__
    for L1b and L2 FITS files and
    `here <https://noaa-goes16.s3.amazonaws.com/index.html>`__ for L1b netCDF files (publicly 
    available Amazon S3 bucket, look on page 2 after all the ABI data).
"""
from .suvi import *  # NOQA
