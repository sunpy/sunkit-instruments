import numpy as np
import pandas as pd
from scipy import interpolate

from astropy import units as u
from astropy.io import fits
from astropy.time import Time

from sunpy import timeseries as ts
from sunpy.data import manager
from sunpy.time import parse_time
from sunpy.util.exceptions import warn_user

__all__ = ["calculate_temperature_em"]

MAX_SUPPORTED_SATELLITE = 17


def calculate_temperature_em(goes_ts, abundance="coronal"):
    """
    This function calculates the isothermal temperature and
    corresponding volume emission measure of the solar soft X-ray
    emitting plasma observed by GOES/XRS.

    These are calculated based on methods described in White et al. 2005 [1]_ (see notes)
    for which the GOES fluxes and channel ratios are used together with look-up tables
    of CHIANTI atomic models to estimate isothermal temperature and emission measure.
    Technically speaking, the method interpolates on the channel flux ratio using
    pre-calcuated tables for the fluxes at a series of temperatures for fixed emission
    measure.

    The method here is almost an exact replica of what is available in SSWIDL,
    namely, goes_chianti_tem.pro, and it has been tested against that for consistency.

    It now also works for the GOES-16 and -17 data, and for the re-processed netcdf
    GOES/XRS files for GOES 13-15.

    Also note that this has only been tested on the high resolutions 1s/2s/3s data of GOES/XRS.

    Parameters
    ----------
    goes_ts : `~sunpy.timeseries.XRSTimeSeries`
        The GOES/XRS timeseries containing the data of both the xrsa and xrsb channels (in units of W/m**2).
    abundance: {"coronal", "photospheric"}, optional
        Which abundances to use for the calculation, the default is "coronal".

    Returns
    -------
    `~sunpy.timeseries.GenericTimeSeries`
        Contains the temperature and emission measure calculated from the input ``goes_ts`` time series.

    Example
    -------
    >>> from sunpy import timeseries as ts
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> from sunkit_instruments import goes_xrs
    >>> goes_ts = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES) # doctest: +REMOTE_DATA
    >>> goes_flare = goes_ts.truncate("2011-06-07 06:20", "2011-06-07 07:30") # doctest: +REMOTE_DATA
    >>> goes_temp_emiss = goes_xrs.calculate_temperature_em(goes_flare) # doctest: +REMOTE_DATA

    Notes
    -----
    This function works with both the NOAA-provided netcdf files, for which the data is given in "true" fluxes
    and with the older FITS files provided by the SDAC, for which the data are scaled to be consistent with GOES-7.
    The routine determines from the `sunpy.timeseries.XRSTimeSeries` metadata
    whether the SWPC scaling factors need to be removed (which are present in the FITS data).

    See also: https://hesperia.gsfc.nasa.gov/goes/goes.html#Temperature/Emission%20Measure

    In regards to the re-processed GOES 8-15 data, please refer to the documentation here:

        * https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/GOES_1-15_XRS_Science-Quality_Data_Readme.pdf
          This includes information of the SWPC scaling factors that were applied to prior data (e.g. the FITS data).

    For the GOES-R data please refer to here:

        * https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/docs/GOES-R_XRS_L2_Data_Readme.pdf

    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005,
        Sol. Phys., 227, 231, DOI: 10.1007/s11207-005-2445-z
    """
    if not isinstance(goes_ts, ts.XRSTimeSeries):
        raise TypeError(
            f"Input time series must be a XRSTimeSeries instance, not {type(goes_ts)}"
        )

    if goes_ts.observatory is None:
        raise ValueError(
            "The GOES satellite number was not found in the input time series"
        )

    satellite_number = int(goes_ts.observatory.split("-")[-1])
    if (satellite_number < 1) or (satellite_number > MAX_SUPPORTED_SATELLITE):
        raise ValueError(
            f"GOES satellite number has to be between 1 and 17, {satellite_number} was found."
        )

    allowed_abundances = ["photospheric", "coronal"]
    if abundance not in allowed_abundances:
        raise ValueError(
            f"The abundance can only be `coronal` or `photospheric`, not {abundance}."
        )
    # Check if GOES-R and whether the primary detector values are given
    if satellite_number >= 16:
        if "xrsa_primary_chan" in goes_ts.columns:
            output = _manage_goesr_detectors(
                goes_ts, satellite_number, abundance=abundance
            )
        else:
            warn_user(
                "No information about primary/secondary detectors in XRSTimeSeries, assuming primary for all"
            )
            output = _chianti_temp_emiss(goes_ts, satellite_number, abundance=abundance)

    # Check if the older files are passed, and if so then the scaling factor needs to be removed.
    # The newer netcdf files now return "true" fluxes so this SWPC factor doesn't need to be removed.
    else:
        if ("Origin" in goes_ts.meta.metas[0]) and (
            goes_ts.meta.metas[0].get("Origin") == "SDAC/GSFC"
        ):
            remove_scaling = True
        else:
            remove_scaling = False
        output = _chianti_temp_emiss(
            goes_ts,
            satellite_number,
            abundance=abundance,
            remove_scaling=remove_scaling,
        )

    return output


@manager.require(
    "goes_chianti_response_table",
    [
        "https://sohoftp.nascom.nasa.gov/solarsoft/gen/idl/synoptic/goes/goes_chianti_response_latest.fits"
    ],
    "cb00c05850e3dc3bbd856eb07c1a372758d689d0845ee591d6e2531afeab0382",
)
def _chianti_temp_emiss(
    goes_ts, satellite_number, secondary=0, abundance="coronal", remove_scaling=False
):
    """
    Calculate isothermal temperature and emission measure from GOES XRS observations.

    This uses the latest formulated response tables including the responses for GOES1-17 to interpolate for temperature or emission measure
    from the measured true fluxes, which is read in from the FITS files goes_chianti_response_latest.fits.

    From the ratio of the two channels the temperature is computed from a spline fit from a lookup response table for 101 temperatures and
    then the emission measure is derived from the temperature and the long flux.

    Parameters
    ----------
    goes_ts : `~sunpy.timeseries.XRSTimeSeries`
        The GOES XRS timeseries containing the data of both the xrsa and xrsb channels (in units of W/m**2).
    sat : `int`
        GOES satellite number.
    secondary: `int`, optional
        Values 0, 1, 2, 3 indicate A1+B1, A2+B1, A1+B2, A2+B2 detector combos for GOES-R.
        Defaults to 0.
    abundance: {"coronal", "photospheric"}, optional
        Which abundances to use for the calculation, the default is "coronal".
    remove_scaling: `bool`, optional
        Checks whether to remove the SWPC scaling factors.
        This is only an issue for the older FITS files for GOES 8-15 XRS.
        Default is `False` as the netcdf files have the "true" fluxes.

    Returns
    -------
    `~sunpy.timeseries.GenericTimeSeries`
        Contains the temperature and emission measure calculated from the input ``goes_ts`` time series.
        The two columns are:
            ``temperature`` : The temperature in MK.
            ``emission `measure` : The volume emission measure.
    Notes
    -----
    Requires goes_chianti_resp.fits produced by goes_chianti_response.pro
    This file contains the pregenerated responses for default coronal and photospheric ion abundances using Chianti version 9.0.1.
    url = "https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits"

    The response table within the fits file starts counting the satellite number at 0.
    For example, for GOES 15 - then 14 is passed to the response table. This is dealt within the code,
    the satellite number to be passed to this function should be the actual GOES satellite number.
    """

    longflux = goes_ts.quantity("xrsb").to(u.W / u.m**2)
    shortflux = goes_ts.quantity("xrsa").to(u.W / u.m**2)

    if "xrsb_quality" in goes_ts.columns:
        longflux[goes_ts.to_dataframe()["xrsb_quality"] != 0] = np.nan
        shortflux[goes_ts.to_dataframe()["xrsa_quality"] != 0] = np.nan

    obsdate = parse_time(goes_ts._data.index[0])

    longflux_corrected = longflux
    # For some reason that I can't find documented anywhere other than in the IDL code,
    # the long channel needs to be scaled by this value for GOES-6 before 1983-06-28.
    if obsdate <= Time("1983-06-28") and satellite_number == 6:
        longflux_corrected = longflux * (4.43 / 5.32)

    shortflux_corrected = shortflux
    # Remove the SWPC scaling factors if needed.
    # The SPWC scaling factors of 0.7 and 0.85 for the XRSA and XSRB channels
    # respectively are documented in the NOAA readme file linked in the docstring.
    if remove_scaling and satellite_number >= 8 and satellite_number < 16:
        longflux_corrected = longflux_corrected / 0.7
        shortflux_corrected = shortflux / 0.85

    # Measurements of short channel flux of less than 1e-10 W/m**2 or
    # long channel flux less than 3e-8 W/m**2 are not considered good.
    # Ratio values corresponding to such fluxes are set to 0.003.
    index = np.logical_or(
        shortflux_corrected < u.Quantity(1e-10 * u.W / u.m**2),
        longflux_corrected < u.Quantity(3e-8 * u.W / u.m**2),
    )
    fluxratio = shortflux_corrected / longflux_corrected
    fluxratio.value[index] = u.Quantity(0.003, unit=u.dimensionless_unscaled)

    # Work out detector index to use from the table response based on satellite number
    # The counting in the table starts at 0, and indexed in an odd way for the GOES-R
    # primary/secondary detectors.
    if satellite_number <= 15:
        sat = satellite_number - 1  # counting starts at 0
    else:
        sat = (
            15 + 4 * (satellite_number - 16) + secondary
        )  # to figure out which detector response table to use (see notes)

    resp_file_name = manager.get("goes_chianti_response_table")
    response_table = fits.getdata(resp_file_name, extension=1)
    rcor = response_table.FSHORT_COR / response_table.FLONG_COR  # coronal
    rpho = response_table.FSHORT_PHO / response_table.FLONG_PHO  # photospheric

    table_to_response_em = 10.0 ** (
        49.0 - response_table["ALOG10EM"][sat]
    )  # for some reason in units of 1e49 (which was to stop overflow errors since 10^49 was
    # too big to express as a standard float in IDL.)

    modeltemp = response_table["TEMP_MK"][sat]
    modelratio = rcor[sat] if abundance == "coronal" else rpho[sat]

    # Calculate the temperature and emission measure:

    # get spline fit to model data to get temperatures given the input flux ratio.
    spline = interpolate.splrep(modelratio, modeltemp, s=0)
    temp = interpolate.splev(fluxratio, spline, der=0)

    modelflux = (
        response_table["FLONG_COR"][sat]
        if abundance == "coronal"
        else response_table["FLONG_PHO"][sat]
    )

    modeltemp = response_table["TEMP_MK"][sat]

    spline = interpolate.splrep(modeltemp, modelflux * table_to_response_em, s=0)
    denom = interpolate.splev(temp, spline, der=0)

    emission_measure = longflux_corrected.value / denom

    goes_times = goes_ts._data.index
    df = pd.DataFrame(
        {"temperature": temp, "emission_measure": emission_measure * 1e49},
        index=goes_times,
    )

    units = {"temperature": u.MK, "emission_measure": u.cm ** (-3)}

    header = {"Info": "Estimated temperature and emission measure"}

    # return a new timeseries with temperature and emission measure
    temp_em = ts.TimeSeries(df, header, units)

    return temp_em


def _manage_goesr_detectors(goes_ts, satellite_number, abundance="coronal"):
    """
    This manages which response to use for the GOES primary and secondary detectors used in the
    observations for the GOES-R satellites (i.e. GOES 16 and 17).

    The GOES 16- and 17- have two detectors for each channel to extend the dynamic range of observations,
    namely XRS-A1, XRS-A2, XRS-B1, XRS-B2, for xrsa and xrsb, respectively.
    During large flares, both the primary and secondary detectors may be used and given that each have a different
    response, we need to match the data at individual times with the correct response.

    Note that the primary channel conditions are values of 0,1,2,3 to indicate A1+B1, A2+B1, A1+B2, A2+B2 detector combos for GOES-R.
    Here, we use the `xrsa{b}_primary_chan` columns to figure out which detectors are used for each timestep.
    """

    # These are the conditions for which detector combinations to use.
    secondary_det_conditions = {0: [1, 1], 1: [2, 1], 2: [1, 2], 3: [2, 2]}

    # here we split the timeseries into sections for which primary or secondary detectors are used
    # and then calculate the response for each and then concatenate them back together.
    outputs = []
    for k in secondary_det_conditions:
        dets = secondary_det_conditions[k]
        (second_ind,) = np.where(
            (goes_ts.quantity("xrsa_primary_chan") == dets[0])
            & (goes_ts.quantity("xrsb_primary_chan") == dets[1])
        )

        goes_split = ts.TimeSeries(goes_ts._data.iloc[second_ind], goes_ts.units)

        if len(goes_split._data) > 0:
            output = _chianti_temp_emiss(
                goes_split, satellite_number, abundance=abundance, secondary=int(k)
            )
            outputs.append(output)

    if len(outputs) > 1:
        full_output = outputs[0].concatenate(outputs[1:])
    else:
        full_output = outputs[0]

    return full_output
