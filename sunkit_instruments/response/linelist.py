"""
CHIANTI line lists with contribution functions (GOFNT).
"""

import os
import logging
import pathlib
import warnings

import numpy as np
import xarray as xr

import astropy.constants as const
import astropy.units as u

__all__ = [
    "LineListEmissionModel",
    "chianti_line_list",
    "get_line_list",
    "line_list_cache_path",
    "line_list_from_emission_model",
]

log = logging.getLogger(__name__)


def chianti_line_list(
    temperature: xr.DataArray,
    density: xr.DataArray = None,
    pressure: xr.DataArray = None,
    abundance: str = None,
    wavelength_range=None,
    minimum_abundance: float = None,
    element_list: list = None,
    ion_list: list = None,
) -> xr.Dataset:
    """
    Generate a line list with contribution functions (GOFNT) using ChiantiPy.

    Parameters
    ----------
    temperature : `xarray.DataArray`
        Temperature array in K with a ``logT`` dimension.
    density : `xarray.DataArray`, optional
        Electron density array in cm^-3. Mutually exclusive with ``pressure``.
    pressure : `xarray.DataArray`, optional
        Electron pressure array in K cm^-3. Mutually exclusive with ``density``.
    abundance : `str`, optional
        CHIANTI abundance name, e.g. ``"sun_coronal_2021_chianti"``.
    wavelength_range : `tuple`
        Two-element (min, max) wavelength range in Angstroms.
    minimum_abundance : `float`, optional
        Minimum elemental abundance to keep.
    element_list : `list`, optional
        List of elements, forwarded to `ChiantiPy.core.bunch`.
    ion_list : `list`, optional
        List of ions, forwarded to `ChiantiPy.core.bunch`. Ignored when
        ``minimum_abundance`` is given.

    Returns
    -------
    `xarray.Dataset`
        Line list with GOFNT and per-transition metadata, restricted to
        ``wavelength_range``.

    Notes
    -----
    The ``XUVTOP`` environment variable must point at a local copy of the
    CHIANTI database. ChiantiPy's ``gui`` default is forced to `False` for
    the calling process, so no ``chiantirc`` file is needed on headless
    systems.
    """
    if (density is None) == (pressure is None):
        msg = "Specify exactly one of density or pressure"
        raise ValueError(msg)
    if not isinstance(temperature, xr.DataArray):
        msg = "temperature must be an xarray.DataArray with a logT dimension"
        raise TypeError(msg)
    wavelength_range = _validate_wavelength_range(wavelength_range)
    if ion_list is not None and minimum_abundance is not None:
        log.warning("minimum_abundance is set, the ion_list will be ignored")

    if "XUVTOP" not in os.environ:
        msg = (
            "The XUVTOP environment variable is not set; ChiantiPy cannot locate the CHIANTI database. "
            "Point it at a local copy of the database, e.g. `export XUVTOP=/path/to/chianti/dbase` "
            "(available from https://www.chiantidatabase.org)."
        )
        raise OSError(msg)

    with warnings.catch_warnings():
        # without a chiantirc file, ChiantiPy evaluates os.path.isfile(False)
        # at import, which raises a RuntimeWarning on Python >= 3.14
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            import ChiantiPy.core as ch
            import ChiantiPy.tools.data as chdata
        except ImportError:
            msg = "ChiantiPy is required for this function, install it with `pip install sunkit-instruments[chianti]`"
            raise ImportError(msg) from None

    # never let ChiantiPy pop GUI selection dialogs (the rcfile default is
    # gui=True, which hangs headless batch jobs)
    chdata.Defaults["gui"] = False

    if density is not None:
        temperature_bc, density_bc = xr.broadcast(temperature, density)
        temperature_flat = temperature_bc.data.reshape(-1)
        density_flat = density_bc.data.reshape(-1)
        extra_coord_name = density.dims[0]
        extra_coord = np.log10(density.data)
    else:
        density = pressure / temperature
        temperature_bc = temperature.broadcast_like(density)
        density_flat = density.data.reshape(-1)
        temperature_flat = temperature_bc.data.reshape(-1)
        extra_coord_name = pressure.dims[0]
        extra_coord = pressure.data

    bunch = ch.bunch(
        temperature_flat,
        density_flat,
        wavelength_range,
        em=1.0,
        abundance=abundance,
        allLines=True,
        keepIons=True,
        minAbund=minimum_abundance,
        ionList=ion_list,
        elementList=element_list,
    )

    return _dataset_from_chianti_bunch(
        bunch,
        temperature,
        temperature_bc,
        {extra_coord_name: extra_coord},
        abundance,
        wavelength_range,
    )


def line_list_cache_path(
    output_dir: pathlib.Path,
    abundance: str,
    wavelength_range,
    *,
    density_dependent: bool = False,
) -> pathlib.Path:
    """
    The canonical cache file path used by `get_line_list` for one abundance.
    """
    wavelength_range = _validate_wavelength_range(wavelength_range)
    prefix = "ll_wvl_eDens" if density_dependent else "ll_wvl"
    lower, upper = (_format_wavelength_bound(bound) for bound in wavelength_range)
    return pathlib.Path(output_dir) / f"{prefix}{lower}_{upper}_{abundance}.ncdf"


def get_line_list(
    *,
    output_dir: pathlib.Path,
    abundance: str,
    wavelength_range,
    temperature: xr.DataArray,
    pressure: xr.DataArray = None,
    density: xr.DataArray = None,
    minimum_abundance: float = None,
    line_list_file: pathlib.Path = None,
    compute_if_missing: bool = True,
) -> xr.Dataset:
    """
    Load the cached CHIANTI line list for one abundance, computing and
    caching it when absent.

    The cache file lives in ``output_dir`` and is named from the wavelength
    range and abundance (``ll_wvl{lo}_{hi}_{abundance}.ncdf``, with an
    ``_eDens`` prefix variant for density-dependent runs); the name is
    re-derived per abundance so multiple abundances never share a cache.
    Writes are atomic (tmp + replace) so concurrent job-array tasks racing
    on the same missing cache cannot observe a half-written file.

    Parameters
    ----------
    output_dir : `pathlib.Path`
        Directory holding the cache files.
    abundance : `str`
        CHIANTI abundance name.
    wavelength_range : array-like
        Two-element (min, max) wavelength range in Angstroms.
    temperature : `xarray.DataArray`
        Temperature grid in K with a ``logT`` dimension.
    pressure, density : `xarray.DataArray`, optional
        Electron pressure grid, or density grid for density-dependent runs
        (mutually exclusive; both forwarded to `chianti_line_list`).
    minimum_abundance : `float`, optional
        Minimum elemental abundance to keep.
    line_list_file : `pathlib.Path`, optional
        Explicit cache file to load, bypassing the derived name.
    compute_if_missing : `bool`, optional
        If `False`, raise instead of computing when the cache is absent
        (useful to make job-array tasks fail fast when the preparation
        step was skipped).

    Returns
    -------
    `xarray.Dataset`
        The line list.
    """
    cache_path = (
        pathlib.Path(line_list_file)
        if line_list_file is not None
        else line_list_cache_path(output_dir, abundance, wavelength_range, density_dependent=density is not None)
    )
    if line_list_file is not None or cache_path.exists():
        log.info(f"Loading line list from {cache_path}")
        return _load_dataset(cache_path)

    if not compute_if_missing:
        msg = f"line-list cache {cache_path} does not exist; run the line-list preparation step first"
        raise FileNotFoundError(msg)

    log.info("Calculating line list")
    line_list = chianti_line_list(
        temperature=temperature,
        pressure=pressure,
        density=density,
        abundance=abundance,
        wavelength_range=wavelength_range,
        minimum_abundance=minimum_abundance,
    )

    # write atomically so concurrent job-array tasks racing on the same
    # missing cache cannot observe a half-written file
    tmp_path = cache_path.with_name(f"{cache_path.name}.tmp{os.getpid()}")
    _save_compressed_netcdf(tmp_path, line_list)
    tmp_path.replace(cache_path)
    return line_list


def _dataset_from_chianti_bunch(bunch, temperature, temperature_bc, extra_coords, abundance, wavelength_range):
    import ChiantiPy
    import ChiantiPy.tools.io as chio

    per_transition = {
        "ion_name": bunch.Intensity["ionS"],
        "wavelength": bunch.Intensity["wvl"],
        "lower_level_label": bunch.Intensity["pretty1"],
        "upper_level_label": bunch.Intensity["pretty2"],
        "lower_level_index": bunch.Intensity["lvl1"],
        "upper_level_index": bunch.Intensity["lvl2"],
        "spectroscopic_name": np.array([bunch.IonInstances[ion].Spectroscopic for ion in bunch.Intensity["ionS"]]),
        "atomic_number": np.array([bunch.IonInstances[ion].Z for ion in bunch.Intensity["ionS"]]),
        "observed": bunch.Intensity["obs"] == "Y",
    }
    line_list = xr.Dataset({name: ("trans_index", values) for name, values in per_transition.items()})

    gofnt_values = bunch.Intensity["intensity"].reshape((*temperature_bc.data.shape, 1, -1))
    line_list["gofnt"] = xr.DataArray(
        gofnt_values,
        dims=(*temperature_bc.dims, "abundance", "trans_index"),
        coords={"logT": np.log10(temperature), **extra_coords, "abundance": np.array([abundance])},
    )
    line_list["logT_peak"] = np.log10(temperature[{"logT": line_list.gofnt.argmax(dim="logT")}])
    line_list["full_name"] = (
        line_list.spectroscopic_name.astype(object) + " " + line_list.wavelength.astype(str).astype(object)
    )

    line_list.attrs["Chiantipy"] = ChiantiPy.__version__
    line_list.attrs["Chianti"] = chio.versionRead()
    line_list.gofnt.attrs["units"] = "erg cm3 / (s sr)"

    # cut down linelist to match wavelength range
    in_range = (line_list.wavelength > wavelength_range[0]) & (line_list.wavelength < wavelength_range[1])
    return line_list.isel(trans_index=in_range)


def _save_compressed_netcdf(path: pathlib.Path, dataset: xr.Dataset) -> None:
    encoding = dict.fromkeys(dataset.data_vars, {"zlib": True, "complevel": 5})
    # h5netcdf rather than netCDF4: the netCDF4 C bindings fail writing
    # variable-length string variables under numpy >= 2.5
    dataset.to_netcdf(path, encoding=encoding, mode="w", engine="h5netcdf")


class LineListEmissionModel:
    """
    A line-list dataset as an emission model for temperature responses.

    Wraps a line list (from `chianti_line_list`, `get_line_list` or
    `line_list_from_emission_model`) so it conforms to the
    `~sunkit_instruments.response.abstractions.EmissionModel` protocol and can
    be passed to `~sunkit_instruments.response.get_temperature_response`.

    The response evaluates the channel's wavelength response at each line
    centre (a delta-function approximation: no thermal broadening, accurate
    to the variation of the wavelength response across a line profile) and
    sums ``gofnt`` x response over all transitions.

    Parameters
    ----------
    line_list : `xarray.Dataset`
        Line list with ``wavelength``, ``gofnt`` and a ``logT`` coordinate.
        Grid dimensions other than ``logT`` and ``trans_index`` (pressure,
        density, abundance) must have size one — select a single entry first.
    gofnt_unit : `str`, optional
        Unit of ``gofnt``; by default read from its ``units`` attribute,
        falling back to the ChiantiPy convention ``erg cm3 / (s sr)``.
    """

    def __init__(self, line_list: xr.Dataset, gofnt_unit: str = None):
        for dim in set(line_list.gofnt.dims) - {"logT", "trans_index"}:
            if line_list.sizes[dim] != 1:
                msg = f"line_list must have a single {dim}; select one before building the model"
                raise ValueError(msg)
            line_list = line_list.isel({dim: 0}, drop=True)
        self.line_list = line_list
        self._gofnt_unit = u.Unit(
            gofnt_unit if gofnt_unit is not None else line_list.gofnt.attrs.get("units", "erg cm3 / (s sr)")
        )

    @property
    def temperature(self):
        return 10 ** self.line_list.logT.data * u.K

    def calculate_temperature_response(self, channel):
        """
        Temperature response of ``channel`` for this line list.

        ``channel`` is compatible with
        `~sunkit_instruments.response.abstractions.AbstractChannel`.
        """
        wave_response = channel.wavelength_response()
        line_wavelength = u.Quantity(self.line_list.wavelength.data, u.AA)
        response_at_lines = np.interp(
            line_wavelength.to_value(channel.wavelength.unit),
            channel.wavelength.value,
            wave_response.value,
            left=0.0,
            right=0.0,
        )
        # gofnt is an energy emissivity; the wavelength response is per photon
        photons_per_erg = 1 / (const.h * const.c / line_wavelength).to_value(u.erg)
        weight = xr.DataArray(response_at_lines * photons_per_erg, dims="trans_index")
        response = (self.line_list.gofnt * weight).sum("trans_index")
        unit = self._gofnt_unit * wave_response.unit * u.photon / u.erg
        return u.Quantity(response.data, unit)


def line_list_from_emission_model(
    model,
    *,
    include_ionization_fraction: bool = True,
    abundance_label: str = "model",
) -> xr.Dataset:
    """
    Build a line list from an emission model's per-transition emissivities.

    Converts the line emissivities of a
    `~sunkit_instruments.response.abstractions.LineEmissionModel` (e.g.
    `synthesizAR.atomic.EmissionModel`, backed by fiasco) into the line-list
    dataset consumed by
    `~sunkit_instruments.response.create_response_function`, matching the
    `chianti_line_list` schema and its gofnt convention
    (``erg cm3 s-1 sr-1``, equilibrium ionization fraction and abundance
    included).

    Parameters
    ----------
    model : `~sunkit_instruments.response.abstractions.LineEmissionModel`
        Emission model exposing ``temperature``, ``density``, iteration over
        ions, and ``get_line_emissivity(ion)`` returning per-transition
        photon emissivities that exclude the ionization fraction (the
        synthesizAR convention).
    include_ionization_fraction : `bool`, optional
        Multiply each ion's equilibrium ionization fraction into the gofnt,
        by default `True`. Set to `False` to keep the model's convention of
        deferring it (e.g. for time-dependent ionization).
    abundance_label : `str`, optional
        Label for the singleton ``abundance`` coordinate.

    Returns
    -------
    `xarray.Dataset`
        Line list with ``wavelength``, ``gofnt``, ``atomic_number``,
        ``ion_name``, ``spectroscopic_name`` and ``full_name`` along
        ``trans_index``, sorted by wavelength.
    """
    gofnt_unit = u.Unit("erg cm3 s-1 sr-1")
    wavelengths = []
    gofnts = []
    atomic_numbers = []
    ion_names = []
    spectroscopic_names = []
    for ion in model:
        wavelength, emissivity = model.get_line_emissivity(ion)
        if wavelength.size == 0 or not np.any(wavelength.value):
            # placeholder entry for ions without level-population data
            continue
        if include_ionization_fraction:
            emissivity = emissivity * ion.ionization_fraction[:, np.newaxis, np.newaxis]
        # photon -> energy emissivity, per steradian (the ChiantiPy gofnt
        # convention used throughout this subpackage)
        photon_energy = (const.h * const.c / wavelength).to(u.erg) / u.photon
        gofnt = (emissivity * photon_energy / (4 * np.pi * u.sr)).to_value(gofnt_unit)
        n_transitions = wavelength.size
        wavelengths.append(wavelength.to_value(u.AA))
        gofnts.append(gofnt)
        atomic_numbers.append(np.full(n_transitions, ion.atomic_number))
        ion_names.append(np.full(n_transitions, getattr(ion, "ion_name", "")))
        spectroscopic_names.append(np.full(n_transitions, ion.ion_name_roman))

    if not gofnts:
        msg = "emission model produced no line emissivities"
        raise ValueError(msg)

    wavelengths = np.concatenate(wavelengths)
    # (logT, density, trans_index) -> (logT, density, abundance, trans_index)
    gofnt = np.concatenate(gofnts, axis=-1)[..., np.newaxis, :]
    atomic_numbers = np.concatenate(atomic_numbers)
    ion_names = np.concatenate(ion_names)
    spectroscopic_names = np.concatenate(spectroscopic_names)
    order = np.argsort(wavelengths)

    line_list = xr.Dataset(
        {
            "wavelength": ("trans_index", wavelengths[order]),
            "gofnt": (("logT", "density", "abundance", "trans_index"), gofnt[..., order]),
            "atomic_number": ("trans_index", atomic_numbers[order]),
            "ion_name": ("trans_index", ion_names[order]),
            "spectroscopic_name": ("trans_index", spectroscopic_names[order]),
        },
        coords={
            "logT": np.log10(model.temperature.to_value(u.K)),
            "density": model.density.to_value(u.cm**-3),
            "abundance": np.array([abundance_label]),
        },
    )
    full_name = [
        f"{name} {wavelength:.3f}"
        for name, wavelength in zip(line_list.spectroscopic_name.values, line_list.wavelength.values, strict=True)
    ]
    line_list["full_name"] = ("trans_index", np.array(full_name, dtype=object))
    line_list.gofnt.attrs["units"] = str(gofnt_unit)
    return line_list


def _load_dataset(path: pathlib.Path) -> xr.Dataset:
    # h5netcdf to match _save_compressed_netcdf (it also reads netCDF4-written
    # HDF5 files); the netCDF4 C bindings misbehave under numpy >= 2.5
    with xr.open_dataset(path, engine="h5netcdf") as dataset:
        return dataset.load()


def _validate_wavelength_range(wavelength_range):
    try:
        lower, upper = wavelength_range
    except (TypeError, ValueError):
        msg = "wavelength_range must contain exactly two values"
        raise ValueError(msg) from None
    return lower, upper


def _format_wavelength_bound(bound):
    return f"{float(bound):g}"
