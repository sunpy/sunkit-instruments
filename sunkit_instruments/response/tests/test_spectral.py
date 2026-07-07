"""
Tests for ``create_response_function``: response integrals, Doppler shifts,
xarray-broadcasting flexibility (parameters carrying extra dims such as a
diffraction order axis), and the windowed scatter-add contaminant path.
None of these tests require ChiantiPy or the CHIANTI database.
"""

import numpy as np
import pytest
import xarray as xr

import astropy.units as u

from sunkit_instruments.response.spectral import create_response_function

NORM = 1e-27
ORDERS = xr.DataArray([1, 2], dims="order", coords={"order": [1, 2]})
VDOP = np.array([-200.0, 0.0, 200.0]) * u.km / u.s

# a CI-like spread of lines across a wide band, on a finer logT grid, to
# exercise the windowed scatter-add over many wavelength bins
WIDE_BAND_LINES = {"wavelength": np.linspace(105.0, 195.0, 30), "logT": np.linspace(4.6, 7.6, 15)}


def synthetic_line_list(n_lines=2, wavelength=None, logT=None):
    """
    Minimal deterministic line list with the fields ``create_response_function``
    consumes: iron lines (near 171 A unless ``wavelength`` is given) with
    Gaussian-in-logT contribution functions.
    """
    wavelength = np.linspace(170.6, 171.4, n_lines) if wavelength is None else np.asarray(wavelength, dtype=float)
    n_lines = wavelength.size
    logT = np.array([5.8, 6.0, 6.2]) if logT is None else np.asarray(logT, dtype=float)
    peaks = np.linspace(1.0, 0.5, n_lines)
    gofnt = peaks[np.newaxis, :] * np.exp(-((logT[:, np.newaxis] - 6.0) ** 2) / 0.02) * 1e-25
    return xr.Dataset(
        {
            "wavelength": ("trans_index", wavelength),
            "atomic_number": ("trans_index", np.full(n_lines, 26)),
            "gofnt": (("logT", "trans_index"), gofnt),
            "full_name": ("trans_index", [f"Fake Fe {i} {w:.3f}" for i, w in enumerate(wavelength)]),
        },
        coords={"logT": logT},
    )


class TestCreateResponseFunctionScalar:
    def test_integral_matches_gofnt(self):
        """
        Each per-line response is a normalized Gaussian in wavelength, so
        integrating over wavelength at vdop=0 must recover gofnt/normalization.
        """
        ll = synthetic_line_list(1)
        resp = create_response_function(
            ll, vdop=VDOP, instrumental_width=0.02, wavelength_range=[170.0, 172.0], wavelength_step_mA=2.0, num_lines_keep=1
        )
        dlam = float(resp.wavelength[1] - resp.wavelength[0])
        integral = (resp.response.sel(vdop=0).isel(line=0) * dlam).sum("wavelength")
        expected = ll.gofnt.isel(trans_index=0) / NORM
        np.testing.assert_allclose(integral.values, expected.values, rtol=1e-3)

    def test_peak_follows_doppler_shift(self):
        ll = synthetic_line_list(1)
        resp = create_response_function(
            ll, vdop=VDOP, instrumental_width=0.02, wavelength_range=[170.0, 172.0], wavelength_step_mA=2.0, num_lines_keep=1
        )
        c_kms = 299792.458
        for v in (-200.0, 200.0):
            peak_wvl = float(resp.wavelength[resp.response.sel(vdop=v).isel(line=0, logT=1).argmax("wavelength")])
            expected = float(ll.wavelength[0]) * (1 + v / c_kms)
            assert abs(peak_wvl - expected) < 3e-3

    def test_scalar_velocity_inputs(self):
        ll = synthetic_line_list(1)
        resp = create_response_function(
            ll,
            vdop=0.0 * u.km / u.s,
            nonthermal_velocity=0.0,
            instrumental_width=0.02,
            wavelength_range=[170.0, 172.0],
            wavelength_step_mA=2.0,
            num_lines_keep=1,
        )
        assert resp.response.sizes["vdop"] == 1
        assert resp.response.sizes["nonthermal_velocity"] == 1

    def test_missing_line_list_fields_raises(self):
        ll = synthetic_line_list(1).drop_vars("atomic_number")
        with pytest.raises(ValueError, match="atomic_number"):
            create_response_function(ll, wavelength_range=[170.0, 172.0], num_lines_keep=1)


class TestCreateResponseFunctionOrderDims:
    def test_instrumental_width_order_dim(self):
        ll = synthetic_line_list(2)
        width = xr.DataArray([0.01, 0.05], dims="order", coords={"order": [1, 2]})
        resp = create_response_function(
            ll, vdop=VDOP, instrumental_width=width, wavelength_range=[170.0, 172.0], wavelength_step_mA=2.0, num_lines_keep=2
        )
        assert "order" in resp.response.dims
        assert resp.response.sizes["order"] == 2
        assert not np.allclose(
            resp.response.isel(order=0).values,
            resp.response.isel(order=1).values,
        )

    def test_equal_widths_give_identical_slices(self):
        ll = synthetic_line_list(1)
        width = xr.DataArray([0.02, 0.02], dims="order", coords={"order": [1, 2]})
        resp = create_response_function(
            ll, vdop=VDOP, instrumental_width=width, wavelength_range=[170.0, 172.0], wavelength_step_mA=2.0, num_lines_keep=1
        )
        np.testing.assert_array_equal(
            resp.response.isel(order=0).values,
            resp.response.isel(order=1).values,
        )

    def test_wavelength_range_order_dependent(self):
        ll = synthetic_line_list(2)
        wavelength_range = [171.0 * 2 / ORDERS - 1.0, 171.0 * 2 / ORDERS + 1.0]
        resp = create_response_function(
            ll, vdop=VDOP, instrumental_width=0.02, wavelength_range=wavelength_range, num_wavelength_bins=64, num_lines_keep=2
        )
        assert set(resp.wavelength_grid.dims) == {"wavelength", "order"}
        assert resp.response.sizes["wavelength"] == 64

    def test_wavelength_range_order_dependent_requires_num_bins(self):
        ll = synthetic_line_list(1)
        wavelength_range = [171.0 * 2 / ORDERS - 1.0, 171.0 * 2 / ORDERS + 1.0]
        with pytest.raises(ValueError, match="num_wavelength_bins"):
            create_response_function(
                ll, vdop=VDOP, instrumental_width=0.02, wavelength_range=wavelength_range, num_lines_keep=1
            )

    def test_window_requires_one_dimensional_wavelength_grid(self):
        ll = synthetic_line_list(3)
        wavelength_range = [171.0 * 2 / ORDERS - 1.0, 171.0 * 2 / ORDERS + 1.0]
        with pytest.raises(ValueError, match="one-dimensional"):
            create_response_function(
                ll,
                wavelength_range=wavelength_range,
                num_wavelength_bins=64,
                num_lines_keep=0,
                window_sigma=8.0,
            )

    def test_contam_sum_with_order_dim(self):
        """num_lines_keep=0 sums every line; the order axis must survive."""
        ll = synthetic_line_list(3)
        width = xr.DataArray([0.01, 0.05], dims="order", coords={"order": [1, 2]})
        resp = create_response_function(
            ll, vdop=VDOP, instrumental_width=width, wavelength_range=[170.0, 172.0], wavelength_step_mA=2.0, num_lines_keep=0
        )
        assert "order" in resp.response.dims
        assert resp.response.sizes["line"] == 1

    def test_effective_area_band_order(self):
        ll = synthetic_line_list(2)
        wavelength = np.linspace(168.0, 174.0, 40)
        ea = xr.DataArray(
            np.full((40, 1), 10.0),
            dims=("wavelength", "band"),
            coords={"wavelength": wavelength, "band": [171]},
        )
        ea = ea * xr.DataArray([1.0, 0.5], dims="order", coords={"order": [1, 2]})
        resp = create_response_function(
            ll,
            vdop=VDOP,
            instrumental_width=0.02,
            wavelength_range=[170.0, 172.0],
            wavelength_step_mA=2.0,
            num_lines_keep=2,
            effective_area=ea,
        )
        assert {"band", "order", "line", "logT", "vdop", "wavelength"} <= set(resp.response.dims)
        # the order axis came in through the effective area, so the two
        # slices must differ by exactly its ratio
        np.testing.assert_allclose(
            resp.response.sel(order=2).values,
            0.5 * resp.response.sel(order=1).values,
            rtol=1e-12,
        )


class TestWindowedContaminants:
    """The windowed scatter-add path (``window_sigma``) must reproduce the
    exact full-grid result."""

    def test_window_matches_full_grid_all_contaminants(self):
        """Wide-band case: num_lines_keep=0, every line summed."""
        ll = synthetic_line_list(**WIDE_BAND_LINES)
        kw = {
            "vdop": np.arange(-200.0, 210.0, 40.0),
            "instrumental_width": 0.0,
            "wavelength_range": [100.0, 200.0],
            "num_lines_keep": 0,
        }
        full = create_response_function(ll, **kw)
        win = create_response_function(ll, window_sigma=8.0, **kw)
        peak = float(np.abs(full["response"]).max())
        xr.testing.assert_allclose(full["response"], win["response"], atol=peak * 1e-9)

    def test_window_with_kept_lines(self):
        """num_lines_keep>0: kept lines stay full-grid, contaminants windowed."""
        ll = synthetic_line_list(**WIDE_BAND_LINES)
        kw = {
            "vdop": np.arange(-100.0, 110.0, 50.0),
            "instrumental_width": 0.01,
            "wavelength_range": [100.0, 200.0],
            "num_lines_keep": 2,
        }
        full = create_response_function(ll, **kw)
        win = create_response_function(ll, window_sigma=8.0, **kw)
        peak = float(np.abs(full["response"]).max())
        xr.testing.assert_allclose(full["response"], win["response"], atol=peak * 1e-9)

    def test_window_off_grid_contaminants_are_zero(self):
        ll = synthetic_line_list(1, wavelength=[500.0])
        resp = create_response_function(
            ll,
            wavelength_range=[100.0, 101.0],
            num_lines_keep=0,
            window_sigma=8.0,
        )
        assert resp.response.dtype.kind != "O"
        assert not resp.response.any()
