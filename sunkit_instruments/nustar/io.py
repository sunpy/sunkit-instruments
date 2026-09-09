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
    The counts and channel data component and the livetime header 
    component.
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
    The effective area data component from the file.
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
    The channel bin data mapping (channel number to energy for the 
    corresponding PHA file) and the RMF & photon channel data.
    """

    with fits.open(file) as hdul:
        channel_data = hdul[1].data
        rmf_and_photon_data = hdul[2].data

    return channel_data, rmf_and_photon_data
