"""
Module containing functions that load NuSTAR files.
"""

from astropy.io import fits

__all__ = ["read_pha", "read_arf", "read_rmf"]


def read_pha(file):
    """
    Read a .pha file and extract data ande header information.

    Parameters
    ----------
    file : `str`, `file-like` or `pathlib.Path`
        A .pha file (see `~astropy.fits.io.open` for details).

    Returns
    -------
    The event list and meta information for the observation.
    """
    with fits.open(file) as hdul:
        data = hdul[1].data
        header = hdul[0].header

    return data, header


def read_arf(file):
    """
    Read a .arf file and extract useful information from it.

    Parameters
    ----------
    file :  `str`, `file-like` or `pathlib.Path`
        A .arf file (see `~astropy.fits.io.open` for details ).

    Returns
    -------
    The low and high boundary of energy bins, and the ancillary response [cm^2] (data['specresp']).
    """
    with fits.open(file) as hdul:
        data = hdul[1].data

    return data


def read_rmf(file):
    """
    Read a .rmf file and extract useful information from it.

    Parameters
    ----------
    file :  `str`, `file-like` or `pathlib.Path`
        A .rmf file (see `~astropy.fits.io.open` for details).

    Returns
    -------
    The low and high boundary of energy bins (data['energ_lo'], data['energ_hi']), number of sub-set channels in the energy
        bin (data['n_grp']), starting index of each sub-set of channels (data['f_chan']),
        number of channels in each sub-set (data['n_chan']), redistribution matrix [counts per photon] (data['matrix']).
    """

    with fits.open(file) as hdul:
        data = hdul[2].data

    return data
