"""
CHIANTI line lists with contribution functions (GOFNT).
"""

import os
import logging
import pathlib
import warnings

import numpy as np
import xarray as xr

__all__ = ["chianti_line_list", "get_line_list", "line_list_cache_path"]

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

    with warnings.catch_warnings():
        # without a chiantirc file, ChiantiPy evaluates os.path.isfile(False)
        # at import, which raises a RuntimeWarning on Python >= 3.14
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            import ChiantiPy
        except ImportError:
            msg = "ChiantiPy is required for this function, install it with `pip install sunkit-instruments[chianti]`"
            raise ImportError(msg) from None

    if "XUVTOP" not in os.environ:
        msg = (
            "The XUVTOP environment variable is not set; ChiantiPy cannot locate the CHIANTI database. "
            "Point it at a local copy of the database, e.g. `export XUVTOP=/path/to/chianti/dbase` "
            "(available from https://www.chiantidatabase.org)."
        )
        raise OSError(msg)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        import ChiantiPy.core as ch
        import ChiantiPy.tools.data as chdata
        import ChiantiPy.tools.io as chio

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
        extra_coord_name,
        extra_coord,
        abundance,
        wavelength_range,
        ChiantiPy.__version__,
        chio.versionRead(),
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
    if line_list_file is not None:
        log.info(f"Loading line list from {line_list_file}")
        return _load_dataset(line_list_file)

    cache_path = line_list_cache_path(output_dir, abundance, wavelength_range, density_dependent=density is not None)

    if cache_path.exists():
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


def _dataset_from_chianti_bunch(
    bunch,
    temperature,
    temperature_bc,
    extra_coord_name,
    extra_coord,
    abundance,
    wavelength_range,
    chiantipy_version,
    chianti_version,
):
    observed = xr.DataArray(
        bunch.Intensity["obs"] == "Y",
        dims="trans_index",
    )

    gofnt_values = bunch.Intensity["intensity"]
    gofnt_values = gofnt_values.reshape((*temperature_bc.data.shape, 1, -1))
    gofnt = xr.DataArray(
        data=gofnt_values,
        dims=(*temperature_bc.dims, "abundance", "trans_index"),
        coords={
            "logT": np.log10(temperature),
            extra_coord_name: extra_coord,
            "abundance": np.array([abundance]),
        },
    )

    ion_names = xr.DataArray(
        data=bunch.Intensity["ionS"],
        dims="trans_index",
    )
    line_wavelengths = xr.DataArray(
        data=bunch.Intensity["wvl"],
        dims="trans_index",
    )
    lower_level_label = xr.DataArray(
        data=bunch.Intensity["pretty1"],
        dims="trans_index",
    )
    upper_level_label = xr.DataArray(
        data=bunch.Intensity["pretty2"],
        dims="trans_index",
    )
    lower_level_index = xr.DataArray(
        data=bunch.Intensity["lvl1"],
        dims="trans_index",
    )
    upper_level_index = xr.DataArray(
        data=bunch.Intensity["lvl2"],
        dims="trans_index",
    )
    spectroscopic_name = xr.DataArray(
        data=np.array([bunch.IonInstances[ion].Spectroscopic for ion in bunch.Intensity["ionS"]]),
        dims="trans_index",
    )
    atomic_number = xr.DataArray(
        data=np.array([bunch.IonInstances[ion].Z for ion in bunch.Intensity["ionS"]]),
        dims="trans_index",
    )
    logT_peak = np.log10(temperature[{"logT": gofnt.argmax(dim="logT")}])
    line_list = xr.Dataset(
        {
            "ion_name": ion_names,
            "wavelength": line_wavelengths,
            "gofnt": gofnt,
            "lower_level_label": lower_level_label,
            "upper_level_label": upper_level_label,
            "lower_level_index": lower_level_index,
            "upper_level_index": upper_level_index,
            "logT_peak": logT_peak,
            "spectroscopic_name": spectroscopic_name,
            "atomic_number": atomic_number,
            "observed": observed,
        }
    )

    full_line_name = line_list.spectroscopic_name.astype(object) + " " + line_list.wavelength.astype(str).astype(object)
    line_list["full_name"] = full_line_name

    # cut down linelist to match wavelength range
    in_range = (line_list.wavelength > wavelength_range[0]) * (line_list.wavelength < wavelength_range[1])

    line_list.attrs["Chiantipy"] = chiantipy_version
    line_list.attrs["Chianti"] = chianti_version

    return line_list.isel(trans_index=in_range)


def _save_compressed_netcdf(path: pathlib.Path, dataset: xr.Dataset) -> None:
    encoding = dict.fromkeys(dataset.data_vars, {"zlib": True, "complevel": 5})
    # h5netcdf rather than netCDF4: the netCDF4 C bindings fail writing
    # variable-length string variables under numpy >= 2.5
    dataset.to_netcdf(path, encoding=encoding, mode="w", engine="h5netcdf")


def _load_dataset(path: pathlib.Path) -> xr.Dataset:
    with xr.open_dataset(path) as dataset:
        return dataset.load()


def _validate_wavelength_range(wavelength_range):
    try:
        lower, upper = wavelength_range
    except (TypeError, ValueError):
        msg = "wavelength_range must contain exactly two values"
        raise ValueError(msg) from None
    return lower, upper


def _format_wavelength_bound(bound):
    return f"{float(np.asarray(bound).item()):g}"
