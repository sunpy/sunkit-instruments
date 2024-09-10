0.5.0 (2023-11-17)
==================

Maintenance release, no new features or bugfixes.

Breaking Changes
----------------

- Increased minimum version of ``sunpy`` to 5.0.0
- Increased minimum version of Python to 3.9

0.4.0 (2023-04-04)
==================

Backwards Incompatible Changes
------------------------------

- This removes the older version of `sunkit_instruments.goes_xrs.calculate_temperature_em` that no longer works for the re-processed netcdf files and new GOES-R data.

  This also removes the `sunkit_instruments.goes_xrs.calculate_radiative_loss_rate` and `sunkit_instruments.goes_xrs.calculate_xray_luminosity` functions that also no longer work in their current form.

  The new `sunkit_instruments.goes_xrs.calculate_temperature_em` function now returns a new sunpy.timeseries.GenericTimeSeries with the temperature and emission measure rather than appending columns to the passed XRSTimeSeries. (`#81 <https://github.com/sunpy/sunkit-instruments/pull/81>`__)


Features
--------

- Create new function (`sunkit_instruments.goes_xrs.calculate_temperature_em`) to calculate the temperature and emission measure from the GOES XRS measurements including the new re-processed XRS data and GOES-R data.

  See :ref:`sphx_glr_generated_gallery_calculate_goes_temperature_and_emission_measure.py` for an example. (`#81 <https://github.com/sunpy/sunkit-instruments/pull/81>`__)


0.3.0 (2022-05-24)
==================

Backwards Incompatible Changes
------------------------------

- Dropped Python 3.7 support. (`#53 <https://github.com/sunpy/sunkit-instruments/pull/53>`__)
- Minimum version of ``sunpy`` is now 4.0 (LTS). (`#61 <https://github.com/sunpy/sunkit-instruments/pull/61>`__)


Features
--------

- Added functions for `SUVI <https://www.swpc.noaa.gov/products/goes-solar-ultraviolet-imager-suvi>`__:

  * :func:`sunkit_instruments.suvi.read_suvi` to read FITS or NetCDF SUVI files.
  * :func:`sunkit_instruments.suvi.files_to_map` to read SUVI files and turn them into :class:`sunpy.map.GenericMap`.
  * :func:`sunkit_instruments.suvi.despike_l1b_file` and :func:`sunkit_instruments.suvi.despike_l1b_array` to despike SUVI L1b files.
  * :func:`sunkit_instruments.suvi.get_response` to get the response function for a given SUVI L1b file or wavelength. (`#61 <https://github.com/sunpy/sunkit-instruments/pull/61>`__)


Bug Fixes
---------

- Fermi pointing file names changed from "_p202_v001" to "_p310_v001" upstream. (`#48 <https://github.com/sunpy/sunkit-instruments/pull/48>`__)


0.2.0 (2021-02-13)
==================

Features
--------

- Add :func:`sunkit_instruments.rhessi.imagecube2map` function to extract `sunpy.map.MapSequence` objects from a RHESSI 4D image cube. (`#35 <https://github.com/sunpy/sunkit-instruments/pull/35>`__)


0.1.0 (2020-09-30)
==================

Features
--------

- Creation of the package with all code from ``sunpy.instr``.
