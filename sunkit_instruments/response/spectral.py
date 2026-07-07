"""
Spectral response functions built from CHIANTI line lists.
"""

import numpy as np
import xarray as xr

import astropy.constants as const
import astropy.units as u

__all__ = ["create_response_function"]

_GAUSSIAN_EXPRESSION = "gofnt_scaled * exp(-0.5 * (shift / width)**2) / gaussian_norm / width"


def create_response_function(
    line_list: xr.Dataset,
    instrumental_width: float = 0,
    normalization: float = 1e-27,
    method: str = "linear",
    vdop: u.Quantity | np.ndarray | list = None,
    nonthermal_velocity: u.Quantity | np.ndarray | list = None,
    wavelength_range: np.ndarray | list = None,
    wavelength_step_mA: float = 4.9,
    num_wavelength_bins: int = None,
    effective_area: xr.DataArray = None,
    num_lines_keep: int = 2,
    band=None,
    window_sigma: float = None,
) -> xr.Dataset:
    """
    Computes a response function as a function of velocity and temperature.

    Parameters
    ----------
    line_list : `xarray.Dataset`
        Line list, can be created with `~sunkit_instruments.response.chianti_line_list`.
    instrumental_width : `float`, optional
        Instrumental width sigma in Angstroms, by default 0. May be a
        `xarray.DataArray` carrying extra dims (e.g. a diffraction order
        axis); the output inherits those dims.
    normalization : `float`, optional
        Normalization in the response function, by default 1e-27.
    method : `str`, optional
        Interpolation method for the effective area, by default "linear".
    vdop : array-like, optional
        Doppler axis array in km/s.
    nonthermal_velocity : array-like, optional
        Nonthermal velocity in km/s.
    wavelength_range : array-like, optional
        Two-element (min, max) wavelength range of the output grid in
        Angstroms. Endpoints may be `xarray.DataArray` objects carrying
        extra dims (requires ``num_wavelength_bins``).
    wavelength_step_mA : `float`, optional
        Wavelength bin size in milli-Angstroms, by default 4.9.
    num_wavelength_bins : `int`, optional
        Number of wavelength bins; overrides ``wavelength_step_mA``.
    effective_area : `xarray.DataArray`, optional
        Effective area with a ``wavelength`` dimension, interpolated onto
        the output wavelength grid and multiplied into the response.
    num_lines_keep : `int`, optional
        Number of lines kept individually (in line-list order); the rest
        are summed into a single "remaining lines" entry, by default 2.
    band : `str`, optional
        Bandpass label used for the summed-lines entry. If `None` it is
        taken from a ``band`` coordinate of ``line_list`` when present.
    window_sigma : `float`, optional
        If set (e.g. 8-10), the summed contaminant response is built by
        evaluating each line's Gaussian only within +/- ``window_sigma`` of
        its (vdop-shifted) centre and scatter-adding it into a preallocated
        full-grid accumulator, instead of computing and summing every line
        across the entire wavelength grid. Same result to the chosen sigma
        cutoff, dramatically faster for wide bands. `None` (default) keeps
        the exact full-grid path.

    Returns
    -------
    `xarray.Dataset`
        Response function (data variable ``response``) with temperature,
        velocity and wavelength axes. When the wavelength grid is
        one-dimensional it is attached as the ``wavelength`` coordinate;
        a grid carrying extra dims (e.g. order) is attached as the
        ``wavelength_grid`` coordinate instead.
    """
    if band is None and "band" in line_list.coords:
        band = line_list.band.astype(str)

    try:
        import periodictable as pt
    except ImportError:
        msg = "periodictable is required for this function, install it with `pip install sunkit-instruments[chianti]`"
        raise ImportError(msg) from None

    try:
        import numexpr as ne
    except ImportError:
        msg = "numexpr is required for this function, install it with `pip install sunkit-instruments[chianti]`"
        raise ImportError(msg) from None

    CC_kms = const.c.to(u.km / u.s).value  # Light speed (km/s)
    CC_ms = const.c.to(u.m / u.s).value  # Light speed (m/s)
    KB = const.k_B.to(u.J / u.K).value  # Boltzmann constant (J/K = kg m^2/K/s^2)
    MP = const.m_p.to(u.kg).value  # Proton mass (kg)

    # logic for single line response function
    if "trans_index" not in line_list.dims:
        line_list = line_list.expand_dims("trans_index")
    _validate_line_list(line_list)

    wavelength_grid = _make_wavelength_grid(line_list, wavelength_range, wavelength_step_mA, num_wavelength_bins)

    if window_sigma is not None and wavelength_grid.ndim != 1:
        msg = "window_sigma only supports one-dimensional wavelength grids"
        raise ValueError(msg)

    if vdop is not None:
        vdop = _velocity_axis(vdop, "vdop")
        line_centers = line_list["wavelength"] * (1 + vdop / CC_kms)
    else:
        line_centers = line_list["wavelength"]
    if nonthermal_velocity is not None:
        nonthermal_velocity = _velocity_axis(nonthermal_velocity, "nonthermal_velocity")

    atomic_mass = _atomic_mass_from_atomic_number(line_list.atomic_number, pt.elements, MP)

    responses = []
    num_lines = line_list.sizes["trans_index"]
    summed_lines_included = False
    gaussian_norm = np.sqrt(2 * np.pi)
    # windowed scatter-add state (only used when window_sigma is not None)
    grid_values = np.asarray(wavelength_grid.values)
    num_wavelengths = wavelength_grid.sizes["wavelength"]
    contaminant_accumulator = None
    contaminant_dims = None
    contaminant_wavelength_axis = None
    contaminant_coords = None
    for i in range(num_lines):
        line_center = line_centers.isel(trans_index=i)

        line_atomic_mass = atomic_mass.isel(trans_index=i)
        thermal_velocity = np.sqrt(KB * 10**line_list.logT / line_atomic_mass)
        thermal_line_width = line_center * thermal_velocity / CC_ms  # (AA) sigma

        if nonthermal_velocity is not None:
            doppler_width = np.sqrt(
                thermal_line_width**2
                + instrumental_width**2
                + (line_list["wavelength"] * (nonthermal_velocity / CC_kms)) ** 2
            )  # (AA)
        else:
            doppler_width = np.sqrt(thermal_line_width**2 + instrumental_width**2)  # (AA)

        if i < num_lines_keep:
            # kept line, stored individually (few lines; full grid is cheap here)
            line_response, gofnt_scaled = _evaluate_gaussian_response(
                ne,
                wavelength_grid,
                line_center,
                doppler_width,
                line_list.gofnt.isel(trans_index=i),
                normalization,
                gaussian_norm,
            )
            line_response = xr.DataArray(line_response, dims=gofnt_scaled.dims, coords=gofnt_scaled.coords)
            line_response = line_response.expand_dims("line")

            line = line_list.full_name.isel(trans_index=i)
            line_response = line_response.assign_coords({"line": line.expand_dims("line")})

            line_response = line_response.assign_coords(
                line_wavelength=line_list.isel(trans_index=i).wavelength.expand_dims("line")
            )
            responses.append(line_response)

        elif window_sigma is not None:
            # windowed scatter-add: evaluate this contaminant's Gaussian only
            # within +/- window_sigma of its (vdop-shifted) centre and add it
            # into a preallocated full-grid accumulator, skipping the ~zero
            # wavelength bins that dominate the full-grid cost.
            summed_lines_included = True
            sigma_max = float(doppler_width.max())
            center_values = np.asarray(line_center.values, dtype=float)
            i0 = max(0, int(np.searchsorted(grid_values, center_values.min() - window_sigma * sigma_max)))
            i1 = min(num_wavelengths, int(np.searchsorted(grid_values, center_values.max() + window_sigma * sigma_max)))
            if i0 < i1:
                block, gofnt_scaled = _evaluate_gaussian_response(
                    ne,
                    wavelength_grid.isel(wavelength=slice(i0, i1)),
                    line_center,
                    doppler_width,
                    line_list.gofnt.isel(trans_index=i),
                    normalization,
                    gaussian_norm,
                )
                if contaminant_accumulator is None:
                    contaminant_dims = gofnt_scaled.dims
                    contaminant_wavelength_axis = contaminant_dims.index("wavelength")
                    shape = list(block.shape)
                    shape[contaminant_wavelength_axis] = num_wavelengths
                    contaminant_accumulator = np.zeros(shape, dtype=block.dtype)
                    contaminant_coords = {d: gofnt_scaled.coords[d] for d in gofnt_scaled.coords}
                elif gofnt_scaled.dims != contaminant_dims:
                    block = np.moveaxis(
                        block, [gofnt_scaled.dims.index(d) for d in contaminant_dims], range(len(contaminant_dims))
                    )
                window = [slice(None)] * contaminant_accumulator.ndim
                window[contaminant_wavelength_axis] = slice(i0, i1)
                contaminant_accumulator[tuple(window)] += block
            elif contaminant_accumulator is None:
                gofnt_scaled = _response_template(
                    wavelength_grid,
                    line_center,
                    doppler_width,
                    line_list.gofnt.isel(trans_index=i),
                    normalization,
                )
                contaminant_dims = gofnt_scaled.dims
                contaminant_wavelength_axis = contaminant_dims.index("wavelength")
                contaminant_accumulator = np.zeros(gofnt_scaled.shape)
                contaminant_coords = {d: gofnt_scaled.coords[d] for d in gofnt_scaled.coords}

        else:
            # full-grid contaminant accumulation (original exact path)
            contaminant_response, gofnt_scaled = _evaluate_gaussian_response(
                ne,
                wavelength_grid,
                line_center,
                doppler_width,
                line_list.gofnt.isel(trans_index=i),
                normalization,
                gaussian_norm,
                accumulator=contaminant_response if summed_lines_included else None,
            )
            summed_lines_included = True

    if summed_lines_included:
        if window_sigma is not None:
            contaminant_response = xr.DataArray(
                contaminant_accumulator,
                dims=contaminant_dims,
                coords=contaminant_coords,
            )
        else:
            contaminant_response = xr.DataArray(
                contaminant_response,
                dims=gofnt_scaled.dims,
                coords=gofnt_scaled.coords,
            )
        contaminant_response = contaminant_response.expand_dims("line")
        # overwrite line coord with correct label
        if band is None:
            contaminant_label = xr.DataArray(
                [f"Remaining {line_list.sizes['trans_index'] - num_lines_keep} lines"], dims="line"
            )
        else:
            contaminant_label = band + xr.DataArray(
                [f" remaining {line_list.sizes['trans_index'] - num_lines_keep} lines"], dims="line"
            )
        contaminant_response = contaminant_response.assign_coords(line=contaminant_label)
        # assign the wvl of the first line in the list to summed lines.  If this is the brightest line, it may well
        # represent where peak intensity occurs
        contaminant_response = contaminant_response.assign_coords(
            {"line_wavelength": line_list.isel(trans_index=0).wavelength.expand_dims("line")}
        )

        # if all lines are summed, no need to broadcast against existing line label for concatenation
        if num_lines_keep == 0:
            contaminant_response["line"] = contaminant_label
        else:
            contaminant_response["line"] = contaminant_label.broadcast_like(line)

        if responses:
            responses.append(contaminant_response)
        else:
            responses = contaminant_response

    responses = xr.concat(responses, dim="line", coords="different", compat="equals")

    ds = xr.Dataset()
    ds["response"] = responses

    if vdop is not None:
        ds = ds.assign_coords({"vdop": vdop})
    if nonthermal_velocity is not None:
        ds = ds.assign_coords({"nonthermal_velocity": nonthermal_velocity})
    response_attrs = ds.response.attrs

    if effective_area is not None:
        interp = effective_area.interp(wavelength=wavelength_grid, method=method)
        ds["response"] = ds.response * interp
        ds.response.attrs.update(response_attrs)
        ds.response.attrs["units"] = str(normalization * u.erg * u.cm**5 / u.s / u.sr / u.angstrom)
    else:
        ds.response.attrs["units"] = str(normalization * u.erg * u.cm**3 / u.s / u.sr / u.angstrom)

    # a 1-D grid becomes the wavelength index coordinate; a grid carrying
    # extra dims (e.g. order) cannot share the dimension's name, so it is
    # attached as wavelength_grid instead
    if wavelength_grid.ndim == 1:
        ds = ds.assign_coords(wavelength=("wavelength", wavelength_grid.data, {"units": str(u.AA)}))
    else:
        ds.coords["wavelength_grid"] = wavelength_grid
        ds.coords["wavelength_grid"].attrs["units"] = str(u.AA)

    ds.attrs["normalization"] = normalization

    return ds


def _velocity_axis(values, dim):
    if isinstance(values, u.Quantity):
        values = values.to(u.km / u.s).value
    values = np.atleast_1d(values)
    return xr.DataArray(values, dims=dim, coords={dim: values})


def _make_wavelength_grid(line_list, wavelength_range, wavelength_step_mA, num_wavelength_bins):
    if wavelength_range is None:
        # shortest and longest line wavelengths, with 1 Angstrom of padding
        wavelength_range = [line_list.wavelength.min().data - 1, line_list.wavelength.max().data + 1]

    if num_wavelength_bins:
        if isinstance(wavelength_range[0], xr.DataArray) or isinstance(wavelength_range[1], xr.DataArray):
            # np.linspace cannot operate on xarray objects directly
            # (https://github.com/pydata/xarray/issues/9043)
            lower, upper = xr.broadcast(xr.DataArray(wavelength_range[0]), xr.DataArray(wavelength_range[1]))
            return xr.DataArray(
                np.linspace(lower.values, upper.values, num=num_wavelength_bins),
                dims=("wavelength", *lower.dims),
                coords=lower.coords,
            )
        return xr.DataArray(
            np.linspace(wavelength_range[0], wavelength_range[1], num=num_wavelength_bins), dims="wavelength"
        )

    if isinstance(wavelength_range[0], xr.DataArray) and wavelength_range[0].ndim > 0:
        msg = (
            "wavelength_range endpoints with extra dims (e.g. order) require num_wavelength_bins; "
            "np.arange only supports scalar endpoints"
        )
        raise ValueError(msg)
    return xr.DataArray(
        np.arange(wavelength_range[0], wavelength_range[1] + wavelength_step_mA / 1e3, wavelength_step_mA / 1e3),  # [A]
        dims="wavelength",
    )


def _validate_line_list(line_list):
    if not isinstance(line_list, xr.Dataset):
        msg = "line_list must be an xarray.Dataset"
        raise TypeError(msg)
    missing = [name for name in ("wavelength", "atomic_number", "gofnt", "full_name") if name not in line_list]
    if missing:
        msg = f"line_list is missing required variables: {', '.join(missing)}"
        raise ValueError(msg)
    if "logT" not in line_list.coords:
        msg = "line_list must include a logT coordinate"
        raise ValueError(msg)
    missing_trans_index = [
        name
        for name in ("wavelength", "atomic_number", "gofnt", "full_name")
        if "trans_index" not in line_list[name].dims
    ]
    if missing_trans_index:
        msg = f"line_list variables must include a trans_index dimension: {', '.join(missing_trans_index)}"
        raise ValueError(msg)


def _atomic_mass_from_atomic_number(atomic_number, elements, proton_mass):
    atomic_mass = atomic_number.copy()
    masses = np.array([elements[int(z)].mass for z in atomic_mass.data.reshape(-1)])
    atomic_mass.data = masses.reshape(atomic_mass.shape) * proton_mass
    return atomic_mass


def _broadcast_response_inputs(wavelength_grid, line_center, doppler_width, gofnt, normalization):
    shift = (wavelength_grid - line_center).broadcast_like(gofnt)
    width, shift = xr.broadcast(doppler_width, shift)
    gofnt_scaled = gofnt.broadcast_like(width) / normalization
    return xr.broadcast(gofnt_scaled, width, shift)


def _response_template(wavelength_grid, line_center, doppler_width, gofnt, normalization):
    return _broadcast_response_inputs(wavelength_grid, line_center, doppler_width, gofnt, normalization)[0]


def _evaluate_gaussian_response(
    numexpr,
    wavelength_grid,
    line_center,
    doppler_width,
    gofnt,
    normalization,
    gaussian_norm,
    *,
    accumulator=None,
):
    gofnt_scaled, width, shift = _broadcast_response_inputs(
        wavelength_grid, line_center, doppler_width, gofnt, normalization
    )
    local_dict = {
        "gofnt_scaled": gofnt_scaled.data,
        "shift": shift.data,
        "width": width.data,
        "gaussian_norm": gaussian_norm,
    }
    expression = _GAUSSIAN_EXPRESSION
    if accumulator is not None:
        local_dict["accumulator"] = accumulator
        expression = f"{expression} + accumulator"
    return numexpr.evaluate(expression, local_dict=local_dict), gofnt_scaled
