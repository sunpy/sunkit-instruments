import numpy as np
import pytest
from numpy.testing import assert_almost_equal, assert_array_equal

import astropy.units as u
from astropy.time import Time
from astropy.units.quantity import Quantity
from sunpy import timeseries
from sunpy.time import TimeRange, is_time_equal, parse_time
from sunpy.util.exceptions import SunpyUserWarning

from sunkit_instruments import goes_xrs as goes
from sunkit_instruments.data.test import get_test_filepath

# Tests for the GOES temperature and emission measure calculations
goes15_fits_filepath = get_test_filepath("go1520110607.fits")  # test old FITS files
goes15_filepath_nc = get_test_filepath(
    "sci_gxrs-l2-irrad_g15_d20170910_v0-0-0_truncated.nc"
)  # test re-processed netcdf files
goes16_filepath_nc = get_test_filepath(
    "sci_xrsf-l2-flx1s_g16_d20170910_v2-1-0_truncated.nc"
)  # test the GOES-R data


@pytest.mark.parametrize(
    "goes_files", [goes15_fits_filepath, goes15_filepath_nc, goes16_filepath_nc]
)
@pytest.mark.remote_data
def test_calculate_temperature_emiss(goes_files):
    goeslc = timeseries.TimeSeries(goes_files)
    goes_temp_em = goes.calculate_temperature_emiss(goeslc)
    # check that it returns a timeseries
    assert isinstance(goes_temp_em, timeseries.GenericTimeSeries)
    # check that both temperature and emission measure in the columns
    assert "temperature" and "emission_measure" in goes_temp_em.columns
    # check units
    assert goes_temp_em.units["emission_measure"].to(u.cm**-3) == 1
    assert goes_temp_em.units["emission_measure"].to(u.MK) == 1
    # check time index isnt changed
    assert np.all(goeslc.time == goes_temp_em.time)


@pytest.mark.remote_data
def test_calculate_temperature_emiss_abundances():
    goeslc = timeseries.TimeSeries(goes15_filepath_nc)
    goes_temp_em = goes.calculate_temperature_emiss(goeslc, abundance="photospheric")
    assert isinstance(goes_temp_em, timeseries.GenericTimeSeries)

    with pytest.raises(ValueError):
        goes.calculate_temperature_emiss(goeslc, abundance="hello")


@pytest.mark.remote_data
def test_calculate_temperature_emiss_errs():
    # check when not a XRS timeseries is passed
    with pytest.raises(TypeError):
        goes.calculate_temperature_emiss([])
    lyra_ts = timeseries.TimeSeries(
        get_test_filepath("lyra_20150101-000000_lev3_std_truncated.fits.gz")
    )
    with pytest.raises(TypeError):
        goes.calculate_temperature_emiss(lyra_ts)


@pytest.mark.remote_data
def test_calculate_temperature_emiss_no_primary_detector_columns_GOESR():
    goeslc = timeseries.TimeSeries(goes16_filepath_nc)
    goeslc_removed_col = goeslc.remove_column("xrsa_primary_chan").remove_column(
        "xrsb_primary_chan"
    )

    with pytest.warns(SunpyUserWarning):
        goes.calculate_temperature_emiss(goeslc_removed_col)


# Test the other GOES-XRS functionality
@pytest.mark.remote_data
def test_goes_event_list():
    # Set a time range to search
    trange = TimeRange("2011-06-07 00:00", "2011-06-08 00:00")
    # Test case where GOES class filter is applied
    result = goes.get_goes_event_list(trange, goes_class_filter="M1")
    assert type(result) == list
    assert type(result[0]) == dict
    assert type(result[0]["event_date"]) == str
    assert type(result[0]["goes_location"]) == tuple
    assert isinstance(result[0]["peak_time"], Time)
    assert isinstance(result[0]["start_time"], Time)
    assert isinstance(result[0]["end_time"], Time)
    assert type(result[0]["goes_class"]) == str
    assert type(result[0]["noaa_active_region"]) == np.int64
    assert result[0]["event_date"] == "2011-06-07"
    assert result[0]["goes_location"] == (54, -21)
    # float errror
    assert is_time_equal(result[0]["start_time"], parse_time((2011, 6, 7, 6, 16)))
    assert is_time_equal(result[0]["peak_time"], parse_time((2011, 6, 7, 6, 41)))
    assert is_time_equal(result[0]["end_time"], parse_time((2011, 6, 7, 6, 59)))
    assert result[0]["goes_class"] == "M2.5"
    assert result[0]["noaa_active_region"] == 11226
    # Test case where GOES class filter not applied
    result = goes.get_goes_event_list(trange)
    assert type(result) == list
    assert type(result[0]) == dict
    assert type(result[0]["event_date"]) == str
    assert type(result[0]["goes_location"]) == tuple
    assert isinstance(result[0]["peak_time"], Time)
    assert isinstance(result[0]["start_time"], Time)
    assert isinstance(result[0]["end_time"], Time)
    assert type(result[0]["goes_class"]) == str
    assert type(result[0]["noaa_active_region"]) == np.int64
    assert result[0]["event_date"] == "2011-06-07"
    assert result[0]["goes_location"] == (54, -21)
    assert is_time_equal(result[0]["start_time"], parse_time((2011, 6, 7, 6, 16)))
    assert is_time_equal(result[0]["peak_time"], parse_time((2011, 6, 7, 6, 41)))
    assert is_time_equal(result[0]["end_time"], parse_time((2011, 6, 7, 6, 59)))
    assert result[0]["goes_class"] == "M2.5"
    assert result[0]["noaa_active_region"] == 11226


def test_flux_to_classletter():
    """
    Test converting fluxes into a class letter.
    """
    fluxes = Quantity(10 ** (-np.arange(9, 2.0, -1)), "W/m**2")
    classesletter = ["A", "A", "B", "C", "M", "X", "X"]
    calculated_classesletter = [goes.flux_to_flareclass(f)[0] for f in fluxes]
    calculated_classnumber = [float(goes.flux_to_flareclass(f)[1:]) for f in fluxes]
    assert_array_equal(classesletter, calculated_classesletter)
    assert_array_equal([0.1, 1, 1, 1, 1, 1, 10], calculated_classnumber)
    # now test the Examples
    assert goes.flux_to_flareclass(1e-08 * u.watt / u.m**2) == "A1"
    assert goes.flux_to_flareclass(0.00682 * u.watt / u.m**2) == "X68.2"
    assert goes.flux_to_flareclass(7.8e-09 * u.watt / u.m**2) == "A0.78"
    assert goes.flux_to_flareclass(0.00024 * u.watt / u.m**2) == "X2.4"
    assert goes.flux_to_flareclass(4.7e-06 * u.watt / u.m**2) == "C4.7"
    assert goes.flux_to_flareclass(6.9e-07 * u.watt / u.m**2) == "B6.9"
    assert goes.flux_to_flareclass(2.1e-05 * u.watt / u.m**2) == "M2.1"


def test_class_to_flux():
    classes = ["A3.49", "A0.23", "M1", "X2.3", "M5.8", "C2.3", "B3.45", "X20"]
    results = Quantity(
        [3.49e-8, 2.3e-9, 1e-5, 2.3e-4, 5.8e-5, 2.3e-6, 3.45e-7, 2e-3], "W/m2"
    )
    for c, r in zip(classes, results):
        assert_almost_equal(r.value, goes.flareclass_to_flux(c).value)


def test_joint_class_to_flux():
    classes = ["A3.49", "A0.23", "M1", "X2.3", "M5.8", "C2.3", "B3.45", "X20"]
    for c in classes:
        assert c == goes.flux_to_flareclass(goes.flareclass_to_flux(c))


# TODO add a test to check for raising error
