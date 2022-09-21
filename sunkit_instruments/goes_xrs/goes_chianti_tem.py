from astropy import units as u 
from astropy.io import fits
from sunpy import timeseries as ts 
from scipy import interpolate
import numpy as np 
from sunpy.time import parse_time
from astropy.time import Time

goes_ts = ts.TimeSeries("/Users/laura.hayes/dev_things/goes_tests/sci_xrsf-l2-flx1s_g16_d20170910_v2-1-0.nc").truncate("2017-09-10 12:00", "2017-09-10 20:00")
goes_15 = ts.TimeSeries("/Users/laura.hayes/dev_things/goes_tests/sci_gxrs-l2-irrad_g15_d20170910_v0-0-0.nc").truncate("2017-09-10 12:00", "2017-09-10 20:00")


def calculate_temperature_em(goes_ts, remove_scaling=False):
    """
    This calculates the temperature and emission measure from the GOES/XRS flux ratios

    """
    # make sure its a GOES-XRS timeseries
    if not isinstance(goes_ts, ts.XRSTimeSeries):
        raise TypeError(f"goests must be a XRSTimeSeries object, not {type(goes_ts)}")

    satellite_number = int(goes_ts.observatory.split("-")[-1])


    longflux = goes_ts.quantity("xrsb").to(u.W/u.m**2)
    shortflux = goes_ts.quantity("xrsa").to(u.W/u.m**2)

    if "xrsb_quality" in goes_ts.columns:
        longflux[goes_ts.to_dataframe()["xrsb_quality"]!=0] = np.nan
        shortflux[goes_ts.to_dataframe()["xrsa_quality"]!=0] = np.nan

    obsdate = parse_time(goes_ts.index[0])

    if obsdate <= Time("1983-06-28") and sat==6:
        longflux_corrected = longflux * (4.43/5.32)
    else:
        longflux_corrected = longflux
        
    if remove_scaling and sat>=8 and sat<16:
        longflux_corrected = longflux_corrected / 0.7
        shortflux_corrected = shortflux / 0.85
    else:
        shortflux_corrected = shortflux
        
    output = goes_chianti_temp_em(shortflux_corrected, longflux_corrected, sat=satellite_number)
    return output

def goes_chianti_temp_em(shortflux_corrected, longflux_corrected, sat=15, abundances="coronal"):
    '''
    Satellite number starts counting at 0. 
    For example, for GOES 15 - you pass 14 to the response table.

    Parameters
    ----------
    fluxratio : `~np.array`
        the ratio of the the short to the long flux channels of GOES-XRS.
    sat: `int`
        the GOES satellite number
    abundances: "coronal" or "photospheric"

    Returns
    -------
    temperature : `np.array` - the temperature in MK
    emission measure : `np.array` - the emission measure in 10^49cm^-3.

    Notes
    -----
    url = "https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits"

    '''
    
    index = np.logical_or(
        shortflux_corrected < u.Quantity(1e-10*u.W/u.m**2),
        longflux_corrected < u.Quantity(3e-8*u.W/u.m**2),
    )
    fluxratio = shortflux_corrected / longflux_corrected
    fluxratio.value[index] = u.Quantity(0.003*u.W/u.m**2)    

    sat = sat-1 # counting starts at 0ß

    resp_file_name = 'goes_chianti_response_latest.fits'
    aa = fits.getdata(resp_file_name, extension=1)
    rcor = aa.FSHORT_COR / aa.FLONG_COR
    rpho = aa.FSHORT_PHO / aa.FLONG_PHO
    
    table_to_response_em = 10.0**(49.-aa["ALOG10EM"][sat])

    modeltemp = aa["TEMP_MK"][sat]
    if abundances=="coronal":
        modelratio = rcor[sat]
    else:
        modelratio = rpho[sat]

    spline = interpolate.splrep(modelratio, modeltemp, s=0)
    #fluxratio = bratio.values
    temp = interpolate.splev(fluxratio, spline, der=0)

    modelflux = aa["FLONG_COR"][sat]
    modeltemp = aa["TEMP_MK"][sat]

    spline = interpolate.splrep(modeltemp, modelflux*table_to_response_em, s=0)
    denom = interpolate.splev(temp, spline, der=0)

    emission_measure = longflux_corrected.value / denom

    return {"temperature":temp*u.MK, "emission measure":emission_measure*1e49*u.cm**(-3)} 







