import numpy as np
import pytest

from sunkit_instruments import suvi

# Test files are all remote data.
pytestmark = pytest.mark.remote_data


def test_suvi_despiking_fits(L1B_FITS):
    _, l1b_fits_data, _ = suvi.read_suvi(
        L1B_FITS,
    )
    despiked_l1b_fits_data = l1b_fits_data
    despiked_l1b_fits_data = suvi.despike_l1b_file(L1B_FITS)
    assert not np.array_equal(l1b_fits_data, despiked_l1b_fits_data)


def test_suvi_despiking_nc(L1B_NC):
    _, l1b_nc_data, _ = suvi.read_suvi(L1B_NC)
    despiked_l1b_nc_data = l1b_nc_data
    despiked_l1b_nc_data = suvi.despike_l1b_file(L1B_NC)
    assert not np.array_equal(l1b_nc_data, despiked_l1b_nc_data)


def test_get_response_nc(L1B_NC):
    l1b_nc_response = suvi.get_response(L1B_NC)
    assert l1b_nc_response["wavelength_channel"] == 171


def test_get_response_fits(L1B_FITS):
    l1b_fits_response = suvi.get_response(L1B_FITS)
    assert l1b_fits_response["wavelength_channel"] == 171


def test_get_response_wavelength():
    response_195 = suvi.get_response(195)
    assert response_195["wavelength_channel"] == 195
