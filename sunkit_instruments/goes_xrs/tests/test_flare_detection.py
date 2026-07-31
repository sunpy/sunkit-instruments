import numpy as np
import pandas as pd
import pytest

import astropy.units as u
from astropy.table import QTable
from astropy.time import Time

from sunpy import timeseries as ts

from sunkit_instruments import goes_xrs as goes

BACKGROUND = 1e-7


def _make_xrs_timeseries(flux, quality=None, cadence="1min", start="2020-01-01"):
    n = len(flux)
    idx = pd.date_range(start, periods=n, freq=cadence)
    data = {"xrsa": flux * 0.5, "xrsb": flux}
    if quality is not None:
        data["xrsb_quality"] = quality
    df = pd.DataFrame(data, index=idx)
    units = {"xrsa": u.W / u.m**2, "xrsb": u.W / u.m**2}
    if quality is not None:
        units["xrsb_quality"] = u.dimensionless_unscaled
    return ts.TimeSeries(df, {}, units, source="xrs")


def _add_flare(flux, start, rise_len, peak_val, decay_tau, decay_len):
    flux = flux.copy()
    flux[start : start + rise_len] = np.linspace(flux[start], peak_val, rise_len)
    decay_idx = np.arange(1, decay_len)
    decay = BACKGROUND + (peak_val - BACKGROUND) * np.exp(-decay_idx / decay_tau)
    flux[start + rise_len : start + rise_len + decay_len - 1] = decay
    return flux


def _make_times(n, cadence="1min", start="2020-01-01"):
    return Time(pd.date_range(start, periods=n, freq=cadence).values, format="datetime64")


def test_find_goes_flares_single_flare():
    flux = np.full(60, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-5, 5, 20)
    goes_ts = _make_xrs_timeseries(flux)

    result = goes.find_goes_flares(goes_ts)

    assert isinstance(result, QTable)
    assert len(result) == 1
    assert isinstance(result["start_time"], Time)
    assert isinstance(result["peak_time"], Time)
    assert isinstance(result["end_time"], Time)
    assert result["goes_class"][0] == "M5"
    assert u.allclose(result["peak_flux"][0], 5e-5 * u.W / u.m**2)
    assert result["start_time"][0] < result["peak_time"][0] < result["end_time"][0]


def test_find_goes_flares_multiple_and_min_class_filter():
    flux = np.full(120, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-5, 5, 20)  # M5
    flux = _add_flare(flux, 50, 5, 5e-7, 3, 10)  # B5, below default C1 threshold
    flux = _add_flare(flux, 80, 8, 2e-4, 6, 25)  # X2
    goes_ts = _make_xrs_timeseries(flux)

    result_default = goes.find_goes_flares(goes_ts)
    assert list(result_default["goes_class"]) == ["M5", "X2"]

    result_all = goes.find_goes_flares(goes_ts, min_class="A1")
    assert list(result_all["goes_class"]) == ["M5", "B5", "X2"]

    result_x_only = goes.find_goes_flares(goes_ts, min_class="X1")
    assert list(result_x_only["goes_class"]) == ["X2"]


def test_find_goes_flares_no_flares():
    flux = np.full(60, BACKGROUND)
    goes_ts = _make_xrs_timeseries(flux)

    result = goes.find_goes_flares(goes_ts)

    assert isinstance(result, QTable)
    assert len(result) == 0
    assert set(result.colnames) == {
        "start_time",
        "peak_time",
        "end_time",
        "goes_class",
        "peak_flux",
    }


def test_find_goes_flares_ignores_bad_quality_data():
    flux = np.full(60, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-5, 5, 20)
    quality = np.zeros(60, dtype=int)
    quality[15] = 1
    goes_ts = _make_xrs_timeseries(flux, quality=quality)

    result = goes.find_goes_flares(goes_ts)

    assert len(result) == 1
    assert result["goes_class"][0] == "M5"


def test_find_goes_flares_wrong_type():
    with pytest.raises(TypeError, match="XRSTimeSeries"):
        goes.find_goes_flares(pd.DataFrame({"xrsb": [1, 2, 3]}))


def test_find_goes_flares_missing_column():
    idx = pd.date_range("2020-01-01", periods=5, freq="1min")
    df = pd.DataFrame({"foo": np.ones(5)}, index=idx)
    goes_ts = ts.TimeSeries(df, {}, {"foo": u.W / u.m**2}, source="xrs")

    with pytest.raises(ValueError, match="xrsb"):
        goes.find_goes_flares(goes_ts)


def test_find_flares_accepts_quantity_flux():
    flux = np.full(60, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-5, 5, 20)
    times = _make_times(60)

    result = goes.find_flares(times, flux * u.W / u.m**2)

    assert len(result) == 1
    assert u.allclose(result["peak_flux"][0], 5e-5 * u.W / u.m**2)


def test_find_flares_requires_quantity_flux():
    flux = np.full(60, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-5, 5, 20)
    times = _make_times(60)

    with pytest.raises(TypeError, match="flux must be an astropy Quantity"):
        goes.find_flares(times, flux)


def test_find_flares_requires_quantity_min_flux_and_flux_threshold():
    flux = (np.full(60, BACKGROUND) * u.W / u.m**2)
    times = _make_times(60)

    with pytest.raises(TypeError, match="min_flux must be an astropy Quantity"):
        goes.find_flares(times, flux, min_flux=1e-6)

    with pytest.raises(TypeError, match="flux_threshold must be an astropy Quantity"):
        goes.find_flares(times, flux, flux_threshold=1e-7)


def test_find_flares_accepts_any_convertible_flux_unit():
    # Same physical values as a typical flare, but in CGS units (erg/s/cm**2)
    # instead of W/m**2.
    flux = np.full(60, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-5, 5, 20)
    times = _make_times(60)

    flux_cgs = (flux * u.W / u.m**2).to(u.erg / u.s / u.cm**2)
    result = goes.find_flares(
        times,
        flux_cgs,
        min_flux=1e-3 * u.erg / u.s / u.cm**2,
        flux_threshold=1e-4 * u.erg / u.s / u.cm**2,
    )

    assert len(result) == 1
    assert result["peak_flux"].unit == u.erg / u.s / u.cm**2
    assert u.allclose(result["peak_flux"][0], (5e-5 * u.W / u.m**2).to(u.erg / u.s / u.cm**2))


def test_find_flares_incompatible_threshold_unit_raises():
    flux = np.full(60, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-5, 5, 20)
    times = _make_times(60)

    with pytest.raises(u.UnitConversionError):
        goes.find_flares(times, flux * u.W / u.m**2, min_flux=1 * u.K)


def test_find_flares_mismatched_lengths():
    times = _make_times(60)
    flux = np.full(30, BACKGROUND)

    with pytest.raises(ValueError, match="same length"):
        goes.find_flares(times, flux)


def test_find_flares_exposes_algorithm_parameters():
    # A short, fast flare: only 3 minutes of rise (2 increasing steps) and 3 of
    # decay (2 decreasing steps).
    flux = np.full(30, BACKGROUND)
    flux[10:13] = np.linspace(BACKGROUND, 5e-5, 3)
    flux[13:18] = BACKGROUND + (5e-5 - BACKGROUND) * np.exp(-np.arange(1, 6) / 1.5)
    times = _make_times(30)

    # With the default min_rise_time/min_decay_time of 4 minutes, this flare is too short to
    # be detected.
    result_default = goes.find_flares(times, flux * u.W / u.m**2)
    assert len(result_default) == 0

    # Lowering min_rise_time/min_decay_time to 2 minutes allows it to be picked up.
    result_relaxed = goes.find_flares(
        times, flux * u.W / u.m**2, min_rise_time=2 * u.min, min_decay_time=2 * u.min
    )
    assert len(result_relaxed) == 1
    assert u.allclose(result_relaxed["peak_flux"][0], 5e-5 * u.W / u.m**2)

    # A larger decay_fraction requires the flux to fall further back towards
    # background before the flare is considered over, pushing the end time later.
    result_deep_decay = goes.find_flares(
        times,
        flux * u.W / u.m**2,
        min_rise_time=2 * u.min,
        min_decay_time=2 * u.min,
        decay_fraction=0.9,
    )
    assert result_deep_decay["end_time"][0] > result_relaxed["end_time"][0]


def test_find_flares_min_flux_filter():
    flux = np.full(60, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-7, 3, 20)  # B5
    times = _make_times(60)
    flux = flux * u.W / u.m**2

    assert len(goes.find_flares(times, flux, min_flux=1e-6 * u.W / u.m**2)) == 0
    assert len(goes.find_flares(times, flux, min_flux=1e-7 * u.W / u.m**2)) == 1
    assert len(goes.find_flares(times, flux, min_flux=None)) == 1


def test_find_flares_min_rise_must_be_at_least_two_cadences():
    times = _make_times(30)
    flux = np.full(30, BACKGROUND) * u.W / u.m**2

    with pytest.raises(ValueError, match="at least two cadences"):
        goes.find_flares(times, flux, min_rise_time=30 * u.s, cadence=1 * u.min)


def test_find_flares_tracks_peak_through_a_short_dip():
    # A small bump with only a 3-minute dip (shorter than the default 4-minute
    # min_decay_time) before a much stronger re-brightening. The dip is too
    # short to lock in the small bump as the flare's peak, so the bigger
    # re-brightening is correctly picked up as the same flare's peak.
    times = _make_times(60)
    flux = np.full(60, BACKGROUND)
    flux[10:15] = np.linspace(BACKGROUND, 2e-6, 5)
    flux[15:19] = np.linspace(2e-6, 1.5e-6, 4)
    flux[19:29] = np.linspace(1.5e-6, 5e-5, 10)
    flux[29:45] = BACKGROUND + (5e-5 - BACKGROUND) * np.exp(-np.arange(1, 17) / 5)

    result = goes.find_flares(times, flux * u.W / u.m**2, min_flux=None)

    assert len(result) == 1
    assert u.allclose(result["peak_flux"][0], 5e-5 * u.W / u.m**2)
    assert result["start_time"][0] == times[10]
    assert result["peak_time"][0] == times[28]


def test_find_flares_interrupted_decay_ends_at_next_flares_start():
    # A flare whose decay is interrupted by a second, much stronger flare
    # before it reaches decay_fraction of the way back to background. The
    # first flare's end_time should be filled in as the second flare's
    # start_time, rather than some arbitrary or missing value, matching the
    # convention used in the official GOES event lists.
    times = _make_times(80)
    flux = np.full(80, BACKGROUND)
    flux[10:20] = np.linspace(BACKGROUND, 2e-5, 10)
    flux[20:25] = np.linspace(2e-5, 1.5e-5, 5)
    flux[25:35] = np.linspace(1.6e-5, 8e-5, 10)
    flux[35:60] = BACKGROUND + (8e-5 - BACKGROUND) * np.exp(-np.arange(1, 26) / 5)

    result = goes.find_flares(times, flux * u.W / u.m**2, min_flux=None)

    assert len(result) == 2
    assert u.allclose(result["peak_flux"][0], 2e-5 * u.W / u.m**2)
    assert u.allclose(result["peak_flux"][1], 8e-5 * u.W / u.m**2)
    assert result["end_time"][0] == result["start_time"][1]


def test_find_flares_naive_accepts_quantity_flux():
    flux = np.full(60, BACKGROUND)
    flux = _add_flare(flux, 10, 10, 5e-5, 5, 20)
    times = _make_times(60)

    result = goes.find_flares_naive(times, flux * u.W / u.m**2)

    assert len(result) == 1
    assert u.allclose(result["peak_flux"][0], 5e-5 * u.W / u.m**2)


def test_find_flares_naive_requires_quantity_flux():
    flux = np.full(60, BACKGROUND)
    times = _make_times(60)

    with pytest.raises(TypeError, match="flux must be an astropy Quantity"):
        goes.find_flares_naive(times, flux)


def test_find_flares_naive_misses_gradual_rise_that_find_flares_catches():
    # A rise whose ratio isn't met within the first min_rise_time (4 samples:
    # 1.18x), but is if the uninterrupted rise is allowed to keep building
    # (1.48x by the 9th sample). find_flares' extending-window search catches
    # this; find_flares_naive, checking only the fixed first window, misses
    # it entirely.
    times = _make_times(60)
    flux = np.full(60, BACKGROUND)
    k = np.arange(16)
    flux[10:26] = BACKGROUND * (1 + 0.06 * k)
    peak = flux[25]
    flux[26:50] = BACKGROUND + (peak - BACKGROUND) * np.exp(-np.arange(1, 25) / 5)
    flux = flux * u.W / u.m**2

    assert len(goes.find_flares(times, flux, min_flux=None)) == 1
    assert len(goes.find_flares_naive(times, flux, min_flux=None)) == 0


def test_find_flares_naive_does_not_patch_interrupted_decay():
    # Same fixture as test_find_flares_interrupted_decay_ends_at_next_flares_start.
    # find_flares patches the first flare's end to the second's start;
    # find_flares_naive instead keeps scanning straight through the second
    # flare's rise and only stops once its own decay happens to cross the
    # first flare's decay level, well after the second flare has started.
    times = _make_times(80)
    flux = np.full(80, BACKGROUND)
    flux[10:20] = np.linspace(BACKGROUND, 2e-5, 10)
    flux[20:25] = np.linspace(2e-5, 1.5e-5, 5)
    flux[25:35] = np.linspace(1.6e-5, 8e-5, 10)
    flux[35:60] = BACKGROUND + (8e-5 - BACKGROUND) * np.exp(-np.arange(1, 26) / 5)
    flux = flux * u.W / u.m**2

    refined = goes.find_flares(times, flux, min_flux=None)
    naive = goes.find_flares_naive(times, flux, min_flux=None)

    assert len(refined) == len(naive) == 2
    assert refined["end_time"][0] == refined["start_time"][1]
    assert naive["end_time"][0] > naive["start_time"][1]


def test_find_flares_rebins_to_requested_cadence():
    # 30 second cadence data, flare rises and decays over a couple of minutes.
    n = 240
    times = _make_times(n, cadence="30s")
    flux = np.full(n, BACKGROUND)
    flux = _add_flare(flux, 20, 20, 5e-5, 10, 40)
    flux = flux * u.W / u.m**2

    result_1min = goes.find_flares(times, flux, cadence=1 * u.min)
    assert len(result_1min) == 1
    assert result_1min["peak_flux"][0] > 1e-5 * u.W / u.m**2

    # Rebinning to a coarser 2-minute cadence still detects the flare, just
    # with fewer, coarser samples backing the same rise/decay durations.
    result_2min = goes.find_flares(
        times, flux, cadence=2 * u.min, min_rise_time=4 * u.min, min_decay_time=4 * u.min
    )
    assert len(result_2min) == 1
    assert result_2min["peak_flux"][0] > 1e-5 * u.W / u.m**2
