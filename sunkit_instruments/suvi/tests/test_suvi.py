import numpy as np
import pytest
import astropy.units as u

from sunkit_instruments import suvi

# Test files are all remote data.
pytestmark = pytest.mark.remote_data


def test_suvi_despiking_fits(L1B_FITS):
    _, l1b_fits_data, _ = suvi.read_suvi(L1B_FITS)
    despiked_l1b_fits_data = suvi.despike_l1b_file(L1B_FITS)
    assert not np.array_equal(l1b_fits_data, despiked_l1b_fits_data)


def test_suvi_despiking_nc(L1B_NC):
    _, l1b_nc_data, _ = suvi.read_suvi(L1B_NC)
    despiked_l1b_nc_data = suvi.despike_l1b_file(L1B_NC)
    assert not np.array_equal(l1b_nc_data, despiked_l1b_nc_data)


def test_get_response_nc(L1B_NC):
    header, _, _ = suvi.read_suvi(L1B_NC)
    ccd_temp_avg = (header["CCD_TMP1"] + header["CCD_TMP2"]) / 2.0 * u.deg_C
    l1b_nc_response = suvi.get_response(L1B_NC)
    assert l1b_nc_response["wavelength_channel"] == 171
    assert l1b_nc_response["ccd_temperature"] == ccd_temp_avg


def test_get_response_fits(L1B_FITS):
    header, _, _ = suvi.read_suvi(L1B_FITS)
    ccd_temp_avg = (header["CCD_TMP1"] + header["CCD_TMP2"]) / 2.0 * u.deg_C
    l1b_fits_response = suvi.get_response(L1B_FITS)
    assert l1b_fits_response["wavelength_channel"] == 171
    assert l1b_fits_response["ccd_temperature"] == ccd_temp_avg


def test_get_response_wavelength():
    response_195 = suvi.get_response(195)
    assert response_195["wavelength_channel"] == 195


def test_get_response_explicit_temperature():
    temp = -70.0 * u.deg_C
    response = suvi.get_response(195, ccd_temperature=temp)
    assert response["ccd_temperature"] == temp


def test_get_response_invalid_temperature():
    temp = -70.0  # Without units
    with pytest.raises(TypeError):
        suvi.get_response(195, ccd_temperature=temp)