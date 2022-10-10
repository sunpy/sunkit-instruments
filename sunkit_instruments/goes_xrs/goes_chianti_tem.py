from astropy import units as u 
from astropy.io import fits
from astropy.time import Time
from sunpy import timeseries as ts 
from sunpy.time import parse_time
from sunpy.data import manager
from scipy import interpolate
import numpy as np 
import pandas as pd
from sunpy.util.exceptions import warn_user

__all__ = ["goes_calculate_temperature_em", "_goes_chianti_temp_em", "_manage_goesr_detectors"]


def goes_calculate_temperature_em(goes_ts, abundance="coronal", remove_scaling=False):
    """
    This calculates the temperature and emission measure from the GOES/XRS flux ratios.

    Parameters
    ----------
    goes_ts : `~sunpy.timeseries.XRSTimeSeries`
        The GOES XRS timeseries containing the data of both the xrsa and xrsb channels (in units of W/m**2).
    abundance: `~str`, {`coronal` | `photospheric`}
        Default 'coronal'. Which abundances to use for the calculation.
    remove_scaling: `Boolean`
        Checks whether to remove the SWPC scaling factors. Default False for "true" fluxes from the netcdf files. 
        This is only an issue and will need to be used for the older FITS files for 8-15 XRS.

    Returns
    -------
    `~sunpy.timeseries.GenericTimeSeries`
        Conatins the temperature and emission measure calculated from the input `goes_ts` TimeSeries.

    Example
    -------
    >>> goes_ts = ts.TimeSeries("sci_xrsf-l2-flx1s_g16_d20170910_v2-1-0.nc").truncate("2017-09-10 12:00", "2017-09-10 20:00")
    >>> output = goes_calculate_temperature_em(goes_ts)

    """
    if not isinstance(goes_ts, ts.XRSTimeSeries):
        raise TypeError(f"goests must be a XRSTimeSeries object, not {type(goes_ts)}")

    satellite_number = int(goes_ts.observatory.split("-")[-1])

    if (satellite_number>=16):  
        if "primary_detector_a" in goes_ts.columns:
            output = _manage_goesr_detectors(goes_ts, satellite_number, abundance=abundance)
        else:
            warn_user("No information about primary/secondary detectors in XRSTimeSeries, assuming primary for all")
            output = _goes_chianti_temp_em(goes_ts, satellite_number, abundance=abundance)
    else: 
        output = _goes_chianti_temp_em(goes_ts, satellite_number, abundance=abundance)
    
    return output


@manager.require('goes_chianti_response_table',
                ['https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits'],
                 '4ca9730fb039e8a04407ae0aa4d5e3e2566b93dfe549157aa7c8fc3aa1e3e04d')

def _goes_chianti_temp_em(goes_ts, satellite_number, secondary=0, abundance="coronal", remove_scaling=False):
    '''


    Parameters
    ----------
    goes_ts : `~sunpy.timeseries.XRSTimeSeries`
        The GOES XRS timeseries containing the data of both the xrsa and xrsb channels (in units of W/m**2).
    sat: `int`
        the GOES satellite number
    abundance: `~str`, {`coronal` | `photospheric`}
        Default 'coronal'. Which abundances to use for the calculation.    abundances: "coronal" or "photospheric"

    Returns
    -------
    `~sunpy.timeseries.GenericTimeSeries`
        Conatins the temperature and emission measure calculated from the input `goes_ts` TimeSeries.
        The two columns are:
            `temperature` : the temperature in MK
            `emission measure` : the emission measure in 10^49cm^-3.

    Notes
    -----
    url = "https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits"

    Satellite number starts counting at 0. 
    For example, for GOES 15 - you pass 14 to the response table.

    '''

    #--------PREP THE DATA----------#
    longflux = goes_ts.quantity("xrsb").to(u.W/u.m**2)
    shortflux = goes_ts.quantity("xrsa").to(u.W/u.m**2)

    if "xrsb_quality" in goes_ts.columns:
        longflux[goes_ts.to_dataframe()["xrsb_quality"]!=0] = np.nan
        shortflux[goes_ts.to_dataframe()["xrsa_quality"]!=0] = np.nan

    obsdate = parse_time(goes_ts._data.index[0])

    if obsdate <= Time("1983-06-28") and satellite_number==6:
        longflux_corrected = longflux * (4.43/5.32)
    else:
        longflux_corrected = longflux
        
    if remove_scaling and satellite_number>=8 and satellite_number<16:
        longflux_corrected = longflux_corrected / 0.7
        shortflux_corrected = shortflux / 0.85
    else:
        shortflux_corrected = shortflux

    index = np.logical_or(
        shortflux_corrected < u.Quantity(1e-10*u.W/u.m**2),
        longflux_corrected < u.Quantity(3e-8*u.W/u.m**2),
    )
    fluxratio = shortflux_corrected / longflux_corrected
    fluxratio.value[index] = u.Quantity(0.003*u.W/u.m**2)    

    #--------WORK OUT DETECTOR INDEX TO USE BASED ON SATELLITE NUMBER----------#
    if satellite_number<=15:
        sat = satellite_number-1 # counting starts at 0
    else:
        sat = 15+4*(satellite_number-16)+secondary # to figure out which detector response table to use (see notes)

    #--------READ RESPONSE TABLE----------#
    resp_file_name = manager.get('goes_chianti_response_table')
    response_table = fits.getdata(resp_file_name, extension=1) 
    rcor = response_table.FSHORT_COR / response_table.FLONG_COR
    rpho = response_table.FSHORT_PHO / response_table.FLONG_PHO
    
    table_to_response_em = 10.0**(49.-response_table["ALOG10EM"][sat])

    modeltemp = response_table["TEMP_MK"][sat]
    if abundance=="coronal":
        modelratio = rcor[sat]
    else:
        modelratio = rpho[sat]

    spline = interpolate.splrep(modelratio, modeltemp, s=0)
    temp = interpolate.splev(fluxratio, spline, der=0)

    modelflux = response_table["FLONG_COR"][sat]
    modeltemp = response_table["TEMP_MK"][sat]

    #--------CALCULATE THE TEMPERATURE AND EMISSION MEASURE----------#
    spline = interpolate.splrep(modeltemp, modelflux*table_to_response_em, s=0)
    denom = interpolate.splev(temp, spline, der=0)

    emission_measure = longflux_corrected.value / denom

    #--------RETURN A SUNPY TIMESERIES----------#
    goes_times = goes_ts._data.index
    df = pd.DataFrame({"temperature":temp, "emission measure":emission_measure*1e49}, 
                       index=goes_times)
    
    units = {"temperature": u.MK, "emission measure": u.cm**(-3)}
    
    header = {"Info":"Estimated temperature and emission measure"}
    
    temp_em = ts.TimeSeries(df, header, units)

    return temp_em


def _manage_goesr_detectors(goes_ts, satellite_number, abundance="coronal"):
    """
    This manages which response to use for the GOES primary detectors used in the
    observations for the GOES-R satellites (i.e. GOES 16 and 17).

    SECONDARY - values 0,1,2,3 indicate A1+B1, A2+B1, A1+B2, A2+B2 detector combos for GOES-R.
    This uses the `primary_detector_a/b` columns to figure out which detectors are used for each timestep.

    """
    
    # these are the conditions of which combinations of detectors to use
    secondary_det_conditions = {0: [1, 1], 1:[2, 1], 2:[1, 2], 3:[2, 2]} 

    outputs = []
    for k in secondary_det_conditions:
        dets = secondary_det_conditions[k]
        second_ind = np.where((goes_ts.quantity("primary_detector_a")==dets[0])\
                             &(goes_ts.quantity("primary_detector_b")==dets[1]))[0]
        
        goes_split = ts.TimeSeries(goes_ts._data.iloc[second_ind], goes_ts.units)

        if len(goes_split._data)>0:
            output = _goes_chianti_temp_em(goes_split, satellite_number, abundance=abundance, secondary=int(k))
            outputs.append(output)

    if len(outputs)>1:
        full_output = outputs[0].concatenate(outputs[1:])
    else:
        full_output = outputs[0]

    return full_output



