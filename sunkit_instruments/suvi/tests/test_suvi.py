import re

import numpy as np
import pytest

from sunkit_instruments import suvi


@pytest.mark.remote_data
def test_suvi_despiking_fits(L1B_FITS):
    _, l1b_fits_data, _ = suvi.read_suvi(
        L1B_FITS,
    )
    despiked_l1b_fits_data = l1b_fits_data
    despiked_l1b_fits_data = suvi.despike_l1b_file(L1B_FITS)
    assert not np.array_equal(l1b_fits_data, despiked_l1b_fits_data)


@pytest.mark.remote_data
def test_suvi_despiking_nc(L1B_NC):
    _, l1b_nc_data, _ = suvi.read_suvi(L1B_NC)
    despiked_l1b_nc_data = l1b_nc_data
    despiked_l1b_nc_data = suvi.despike_l1b_file(L1B_NC)
    assert not np.array_equal(l1b_nc_data, despiked_l1b_nc_data)


@pytest.mark.remote_data
def test_get_response_nc(L1B_NC):
    l1b_nc_response = suvi.get_response(L1B_NC)
    assert l1b_nc_response["wavelength_channel"] == 171


@pytest.mark.remote_data
def test_get_response_fits(L1B_FITS):
    l1b_fits_response = suvi.get_response(L1B_FITS)
    assert l1b_fits_response["wavelength_channel"] == 171


@pytest.mark.remote_data
def test_get_response_wavelength():
    response_195 = suvi.get_response(195)
    assert response_195["wavelength_channel"] == 195


@pytest.mark.parametrize("spacecraft", [16, 17, 18, 19])
def test_get_response_spacecraft_number(spacecraft):
    response_195 = suvi.get_response(195, spacecraft=spacecraft)
    assert response_195["wavelength_channel"] == 195


def test_get_response_bad_spacecraft_number():
    with pytest.raises(ValueError, match=re.escape("Invalid spacecraft: 0 Valid spacecraft are: [16, 17, 18, 19]")):
        suvi.get_response(195, spacecraft=0)
