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
