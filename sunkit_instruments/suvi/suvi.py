from pathlib import Path

import numpy as np
from scipy import interpolate
from scipy.ndimage import gaussian_filter

import sunpy.map
from astropy import units as u

from sunkit_instruments.suvi._variables import (
    FILTER_SETUP,
    FLIGHT_MODEL,
    VALID_SPACECRAFT,
    VALID_WAVELENGTH_CHANNELS,
)

__all__ = [
    "despike_l1b_file",
    "despike_l1b_array",
    "get_response",
]

PATH_TO_FILES = Path(__file__).parent / "data"


def _despike(image, dqf_mask, filter_width):
    """
    Helper function to do the actual despiking.

    Parameters
    ----------
    image : `numpy.ndarray`
        Array to despike.
    dqf_mask : `numpy.ndarray`
        Data quality flags array.
    filter_width: `int`, optional.
        The filter width for the Gaussian filter, default is 7.
        If NaNs are still present in the despiked image, try increasing this value.
    """
    image_with_nans = np.copy(image)
    image_with_nans[np.where(dqf_mask == 4)] = np.nan
    indices = np.where(np.isnan(image_with_nans))
    image_gaussian_filtered = gaussian_filter(image, filter_width)
    despiked_image = np.copy(image_with_nans)
    despiked_image[indices] = image_gaussian_filtered[indices]
    return despiked_image


def despike_l1b_file(filename, filter_width=7):
    """
    Despike SUVI L1b data and return a despiked `~sunpy.map.Map`.

    .. warning::
        The despiking relies on the presence of the data quality
        flags (DQF) in the first extension of a SUVI L1b FITS file.
        Early in the mission, the DQF extension was not present
        yet, so the despiking cannot be done with this function
        for those early files.

    Parameters
    ----------
    filename: `str`
        File to despike.
    filter_width: `int`, optional.
        The filter width for the Gaussian filter, default is 7.
        If NaNs are still present in the despiked image, try increasing this value.

    Returns
    -------
    `~sunpy.map.Map`
        The despiked L1b image as a `~sunpy.map.Map`.
    """
    # Avoid circular import
    from sunkit_instruments.suvi.io import read_suvi

    header, image, dqf_mask = read_suvi(filename)
    despiked_image = _despike(image, dqf_mask, filter_width)
    return sunpy.map.Map(despiked_image, header)


def despike_l1b_array(data, dqf, filter_width=7):
    """
    Despike SUVI L1b data and return a despiked `numpy.ndarray`.

    Parameters
    ----------
    data : `numpy.ndarray`
        Array to despike.
    dqf : `numpy.ndarray`
        Data quality flags array.
    filter_width: `int`, optional.
        The filter width for the Gaussian filter, default is 7.
        If NaNs are still present in the despiked image, try increasing this value.

    Returns
    -------
    `numpy.ndarray`
        The despiked L1b image as a numpy array.
    """
    return _despike(data, dqf, filter_width)


def get_response(request, spacecraft=16, ccd_temperature=-60.0, exposure_type="long"):
    """
    Get the SUVI instrument response for a specific wavelength channel,
    spacecraft, CCD temperature, and exposure type.

    ``request`` can either be an L1b filename (FITS or netCDF), in which case all of those
    parameters are read automatically from the metadata, or the parameters
    can be passed manually, with ``request`` specifying the desired wavelength
    channel.

    Parameters
    ----------
    request: `str` or {94 | 131 | 171 | 195 | 284 | 304}.
        Either an L1b filename (FITS or netCDF), or an integer
        specifying the wavelength channel.
    spacecraft: `int`, optional.
        Which GOES spacecraft, default is 16.
    ccd_temperature: `float`, optional.
        The CCD temperature, in degrees Celsius, default is -60.
        Needed for getting the correct gain number.
    exposure_type: {"long" | "short" | "short_flare"}, optional.
        The exposure type of the SUVI image.
        The exposure type is needed for the correct focal plane
        filter selection.

        Can be:
        * "long", "short", "short_flare" for 94 and 131
        * "long", "short_flare" for 171, 195, 284, and 304.

    Returns
    -------
    `dict`
        The instrument response information.
        Keys:

        * "wavelength"
        * "effective_area"
        * "response"
        * "wavelength_channel"
        * "spacecraft"
        * "ccd_temperature"
        * "exposure_type"
        * "flight_model"
        * "gain"
        * "solid_angle"
        * "geometric_area"
        * "filter_setup"
    """
    # Avoid circular import
    from sunkit_instruments.suvi.io import read_suvi

    if isinstance(request, str):
        header, _, _ = read_suvi(request)
        wavelength_channel = int(header["WAVELNTH"])
        spacecraft = int(header["TELESCOP"].replace(" ", "").replace("G", ""))
        ccd_temperature = (header["CCD_TMP1"] + header["CCD_TMP2"]) / 2.0
        exposure_type = "_".join(
            header["SCI_OBJ"].replace(" ", "").split(sep="_")[3:]
        ).replace("_exposure", "")
    elif isinstance(request, int):
        wavelength_channel = request
    else:
        raise TypeError(
            f"Input not recognized, must be str for filename or int for wavelength channel, not {type(request)}"
        )

    if wavelength_channel not in VALID_WAVELENGTH_CHANNELS:
        raise ValueError(
            f"Invalid wavelength channel: {wavelength_channel}"
            f"Valid wavelength channels are: {VALID_WAVELENGTH_CHANNELS}"
        )
    if spacecraft not in VALID_SPACECRAFT:
        raise ValueError(
            f"Invalid spacecraft: {spacecraft}"
            f"Valid spacecraft are: {VALID_SPACECRAFT}"
        )

    eff_area_file = (
        PATH_TO_FILES
        / f"SUVI_{FLIGHT_MODEL[spacecraft]}_{wavelength_channel}A_eff_area.txt"
    )
    gain_file = PATH_TO_FILES / f"SUVI_{FLIGHT_MODEL[spacecraft]}_gain.txt"

    eff_area = np.loadtxt(eff_area_file, skiprows=12)
    wave = eff_area[:, 0] * u.Angstrom
    if FILTER_SETUP[wavelength_channel][exposure_type]["FW2"] == "open":
        effective_area = eff_area[:, 1] * u.cm * u.cm
    else:
        effective_area = eff_area[:, 2] * u.cm * u.cm

    gain_table = np.loadtxt(gain_file, skiprows=7)
    temp_x = gain_table[:, 0]
    gain_y = gain_table[:, 1]
    gain_vs_temp = interpolate.interp1d(temp_x, gain_y)
    gain = gain_vs_temp(ccd_temperature)

    geometric_area = 19.362316 * u.cm * u.cm
    solid_angle = ((2.5 / 3600.0 * (np.pi / 180.0)) ** 2.0) * u.sr
    master_e_per_phot = (
        (6.626068e-34 * (u.J / u.Hz)) * (2.99792458e8 * (u.m / u.s))
    ) / (wave.to(u.m) * ((u.eV.to(u.J, 3.65)) * u.J))
    response = effective_area * (master_e_per_phot / gain)

    response_info = {
        "wavelength": wave,
        "effective_area": effective_area * (u.ct / u.ph),
        "response": response,
        "wavelength_channel": wavelength_channel,
        "spacecraft": "GOES-" + str(spacecraft),
        "ccd_temperature": ccd_temperature * u.deg_C,
        "exposure_type": exposure_type,
        "flight_model": FLIGHT_MODEL[spacecraft],
        "gain": float(gain),
        "solid_angle": solid_angle,
        "geometric_area": geometric_area,
        "filter_setup": FILTER_SETUP[wavelength_channel][exposure_type],
    }
    return response_info
