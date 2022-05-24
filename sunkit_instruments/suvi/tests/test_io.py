import numpy as np
import pytest

import sunpy.map
from astropy.io import fits
from sunpy.map.sources.suvi import SUVIMap
from sunpy.util.exceptions import SunpyUserWarning

from sunkit_instruments import suvi

# Test files are all remote data.
pytestmark = pytest.mark.remote_data


def test_read_suvi_l1b_nc(L1B_NC):
    l1b_nc_header, l1b_nc_data, l1b_nc_dqf = suvi.read_suvi(L1B_NC)
    assert isinstance(l1b_nc_header, fits.header.Header)
    assert l1b_nc_data.shape == (1280, 1280)
    assert l1b_nc_dqf.shape == (1280, 1280)


def test_read_suvi_l1b_fits(L1B_FITS):
    l1b_fits_header, l1b_fits_data, l1b_fits_dqf = suvi.read_suvi(L1B_FITS)
    assert isinstance(l1b_fits_header, fits.header.Header)
    assert l1b_fits_data.shape == (1280, 1280)
    assert l1b_fits_dqf.shape == (1280, 1280)


def test_read_suvi_l2_composite(L2_COMPOSITE):
    l2_header, l2_data, _ = suvi.read_suvi(L2_COMPOSITE)
    assert isinstance(l2_header, fits.header.Header)
    assert l2_data.shape == (1280, 1280)


def test_suvi_fix_l1b_header(L1B_FITS):
    header = suvi.io._fix_l1b_header(L1B_FITS)
    assert isinstance(header, fits.header.Header)


def test_files_to_map_l1b_nc(L1B_NC):
    one = suvi.files_to_map(L1B_NC)
    collection = suvi.files_to_map([L1B_NC, L1B_NC, L1B_NC, L1B_NC])
    collection_despike = suvi.files_to_map(
        [L1B_NC, L1B_NC, L1B_NC, L1B_NC], despike_l1b=True
    )

    assert isinstance(one, sunpy.map.GenericMap)
    assert isinstance(collection, sunpy.map.MapSequence)
    assert isinstance(collection_despike, sunpy.map.MapSequence)
    np.testing.assert_equal(collection[0].data, one.data)
    assert not np.array_equal(collection[0].data, collection_despike[0].data)

    with pytest.warns(SunpyUserWarning, match="List of data/headers is empty."):
        suvi.files_to_map([L1B_NC, L1B_NC, L1B_NC, L1B_NC], only_short_exposures=True)


def test_files_to_map_l1b_fits(L1B_FITS):
    one = suvi.files_to_map(L1B_FITS)
    collection = suvi.files_to_map([L1B_FITS, L1B_FITS, L1B_FITS, L1B_FITS])
    collection_despike = suvi.files_to_map(
        [L1B_FITS, L1B_FITS, L1B_FITS, L1B_FITS], despike_l1b=True
    )

    assert isinstance(one, sunpy.map.GenericMap)
    assert isinstance(collection, sunpy.map.MapSequence)
    assert isinstance(collection_despike, sunpy.map.MapSequence)
    np.testing.assert_equal(collection[0].data, one.data)
    assert not np.array_equal(collection[0].data, collection_despike[0].data)

    with pytest.warns(SunpyUserWarning, match="List of data/headers is empty."):
        suvi.files_to_map(
            [L1B_FITS, L1B_FITS, L1B_FITS, L1B_FITS], only_short_exposures=True
        )


def test_files_to_map_nc(L1B_NC):
    l1b_nc_map = suvi.files_to_map(L1B_NC)
    assert isinstance(l1b_nc_map, SUVIMap)


def test_files_to_map_fit(L1B_FITS):
    l1b_fits_map = suvi.files_to_map(L1B_FITS)
    assert isinstance(l1b_fits_map, SUVIMap)


def test_files_to_map_l2_composite(L2_COMPOSITE):
    l2_map = suvi.files_to_map(L2_COMPOSITE)
    assert isinstance(l2_map, SUVIMap)
