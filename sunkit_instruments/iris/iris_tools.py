
"""Some IRIS instrument tools."""

from __future__ import division

import datetime
# import warnings
import os.path

import numpy as np
# import astropy
import astropy.units as u
from astropy.units.quantity import Quantity
# from astropy.modeling import fitting
from astropy.modeling.models import custom_model
from astropy import constants
# from astropy.time import Time
from astropy.table import Table
import scipy
import scipy.io
# from scipy import ndimage
from scipy import interpolate
from sunpy.time import parse_time
import sunpy.util.config
from sunpy.util.net import check_download_file
# from ndcube import NDCube


# Define some properties of IRIS detectors.  Source: IRIS instrument paper.
DETECTOR_GAIN = {"NUV": 18., "FUV": 6., "SJI": 18.}
DETECTOR_YIELD = {"NUV": 1., "FUV": 1.5, "SJI": 1.}
SJI_DEFAULT_BSCALE = 0.25
SJI_DEFAULT_BZERO = 7992.0
DN_UNIT = {
    "NUV": u.def_unit("DN_IRIS_NUV",
                      DETECTOR_GAIN["NUV"] / DETECTOR_YIELD["NUV"]*u.photon),
    "FUV": u.def_unit("DN_IRIS_FUV",
                      DETECTOR_GAIN["FUV"]/DETECTOR_YIELD["FUV"]*u.photon),
    "SJI": u.def_unit("DN_IRIS_SJI",
                      DETECTOR_GAIN["SJI"]/DETECTOR_YIELD["SJI"]*u.photon),
    "SJI_UNSCALED": u.def_unit("DN_IRIS_SJI_UNSCALED", u.ct)}
# Define an equivalency between SJI and SJI_UNSCALED units
SJI_SCALING = [(DN_UNIT["SJI"],
                DN_UNIT["SJI_UNSCALED"],
                lambda x: (x - SJI_DEFAULT_BZERO) / SJI_DEFAULT_BSCALE,
                lambda x: x * SJI_DEFAULT_BSCALE + SJI_DEFAULT_BZERO)]

READOUT_NOISE = {"NUV": 1.2*DN_UNIT["NUV"],
                 "FUV": 3.1*DN_UNIT["FUV"],
                 "SJI": 1.2*DN_UNIT["SJI"]}
RADIANCE_UNIT = u.erg / u.cm ** 2 / u.s / u.steradian / u.Angstrom
SLIT_WIDTH = 0.33*u.arcsec

IRIS_RESPONSE_REMOTE_PATH = "https://sohowww.nascom.nasa.gov/solarsoft/iris/response/"
RESPONSE_VERSION_FILENAMES = {"1": "iris_sra_20130211.geny",
                              "2": "iris_sra_20130715.geny",
                              "3": "iris_sra_c_20150331.geny",
                              "4": "iris_sra_c_20161022.geny"}

# Define some custom error messages.
APPLY_EXPOSURE_TIME_ERROR = ("Exposure time correction has probably already "
                             "been applied since the unit already includes "
                             "inverse time. To apply exposure time correction "
                             "anyway, set 'force' kwarg to True.")
UNDO_EXPOSURE_TIME_ERROR = ("Exposure time correction has probably already "
                            "been undone since the unit does not include "
                            "inverse time. To undo exposure time correction "
                            "anyway, set 'force' kwarg to True.")

# Define whether IRIS WCS is 0 or 1 origin based.
WCS_ORIGIN = 1


def get_iris_response(time_obs=None, pre_launch=False, response_file=None, response_version=None,
                      force_download=False, *args, **kwargs):
    """Returns IRIS response structure.

    One and only one of pre_launch, response_file and response_version must be set.

    Parameters
    ----------
    time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2
        Observation times of the datapoints.
        Must be in the format of, e.g.,
        time_obs=parse_time('2013-09-03', format='utime'),
        which yields 1094169600.0 seconds in value.
        The argument time_obs is ignored for versions 1 and 2.

    pre_launch: `bool`
        Equivalent to setting response_version=2.  Cannot be set
        simultaneously with response_file kwarg. Default=False.

    response_file: `int`
        Version number of effective area file to be used.  Cannot be set
        simultaneously with pre_launch kwarg.  Default=latest.

    response_version : `int`
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=4.

    Returns
    -------
    iris_response: `dict`
        Various parameters regarding IRIS response or effective area structure.
        Includes the following keys:
        date_obs: `astropy.time.Time`
        lambda: `astropy.units.Quantity`
        area_sg: `astropy.units.Quantity`
        name_sg: `str`
        dn2phot_sg: `tuple` of length 2
        area_sji: `astropy.units.Quantity`
        name_sji: `str`
        dn2phot_sji:  `tuple` of length 4
        comment: `str`
        version: `int`
        version_date: `astropy.time.Time`

    Notes
    -----
    This routine does calculate time dependent effective areas using
    version 3 and above of the response functions as is done in the SSW version
    of this code. This code has been updated to calculate time-dependent
    effective areas.

    """
    # Ensures the file exits in the path given.
    if response_file is not None:
        if not(os.path.isfile(response_file)):
            raise KeyError("Not a valid file path")

    # Ensure conflicting kwargs are not set.
    if pre_launch:
        pre_launch_set = True
    else:
        pre_launch_set = False
    if response_file:
        response_file_set = True
    else:
        response_file_set = False
    if response_version:
        response_version_set = True
    else:
        response_version_set = False
    if pre_launch_set + response_file_set + response_version_set > 1:
        raise ValueError("One and only one of kwargs pre_launch, response_file "
                         "and response_version must be set.")
    # If pre_launch set, define response_version to 2.
    if pre_launch:
        response_version = 2
    if (response_file_set is False) and (response_version_set is False):
        response_version = 4
    # If response_file not set, define appropriate response file
    # based on version.
    if not response_file:
        try:
            response_filename = RESPONSE_VERSION_FILENAMES[str(response_version)]
        except KeyError:
            raise KeyError("Version number not recognized.")
        # Define the directory in which the response file should exist
        # to be the sunpy download directory.
        config = sunpy.util.config.load_config()
        download_dir = config.get('downloads', 'download_dir')
        # Check response file exists in download_dir.  If not, download it.
        check_download_file(response_filename, IRIS_RESPONSE_REMOTE_PATH, download_dir, replace=force_download)
        # Define response file as path + filename.
        response_file = os.path.join(download_dir, response_filename)
    # Read response file and store in a dictionary.
    raw_response_data = scipy.io.readsav(response_file)
    iris_response = dict([(name, raw_response_data["p0"][name][0]) for name in raw_response_data["p0"].dtype.names])
    # Convert some properties to more convenient types.
    iris_response["LAMBDA"] = Quantity(iris_response["LAMBDA"], unit=u.nm)
    iris_response["AREA_SG"] = Quantity(iris_response["AREA_SG"], unit=u.cm**2)
    iris_response["AREA_SJI"] = Quantity(iris_response["AREA_SJI"], unit=u.cm**2)
    iris_response["GEOM_AREA"] = Quantity(iris_response["GEOM_AREA"], unit=u.cm**2)
    iris_response["VERSION"] = iris_response["VERSION"]
    # Convert some properties not found in version below version 3 to
    # more convenient types.
    if int(iris_response["VERSION"]) > 2:
        # If DATE_OBS has a value, convert to `astropy.time.Time`, else set to None.
        try:
            iris_response["DATE_OBS"] = parse_time(iris_response["DATE_OBS"], format = 'utime')
        except Exception:
            iris_response["DATE_OBS"] = None

        time_obs = np.array([time_obs.value])

        # Convert C_F_TIME to array of time objects while conserving shape.
        c_f_time = np.empty(iris_response["C_F_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_F_TIME"]):
            for j, t in enumerate(row):
                c_f_time[i][j] = parse_time(t, format = 'utime')
        iris_response["C_F_TIME"] = c_f_time

        # Convert C_F_LAMBDA to Quantity.
        iris_response["C_F_LAMBDA"] = Quantity(iris_response["C_F_LAMBDA"], unit="nm")

        # Convert C_N_TIME to array of time objects while
        # conserving shape.
        c_n_time = np.empty(iris_response["C_N_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_N_TIME"]):
            for j, t in enumerate(row):
                c_n_time[i][j] = parse_time(t, format = 'utime')
        iris_response["C_N_TIME"] = c_n_time

        # Convert C_N_LAMBDA to Quantity.
        iris_response["C_N_LAMBDA"] = Quantity(iris_response["C_N_LAMBDA"], unit="nm")

        # Convert C_S_TIME to array of time objects while
        # conserving shape.
        c_s_time = np.empty(iris_response["C_S_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_S_TIME"]):
            for j, column in enumerate(row):
                for k, t in enumerate(column):
                    c_s_time[i][j][k] = parse_time(t, format = 'utime')
        iris_response["C_S_TIME"] = c_s_time

        # Convert DATE in ELEMENTS array to array of time objects.
        for i, t in enumerate(iris_response["ELEMENTS"]["DATE"]):
            iris_response["ELEMENTS"]["DATE"][i] = parse_time(t, format = 'iso')
            # Convert VERSION_DATE to time object.
            iris_response["VERSION_DATE"] = parse_time(iris_response["VERSION_DATE"])

    if int(iris_response["VERSION"]) < 2:
        # Change DATE tag in data with version < 2 to VERSION_DATE to
        # be consistent with more recent versions.
        iris_response["VERSION_DATE"] = Time(datetime.datetime(int(iris_response["DATE"][0:4]),
                int(iris_response["DATE"][4:6]),
                int(iris_response["DATE"][6:8])))
        del(iris_response["DATE"])

    if int(iris_response["VERSION"]) > 2 and time_obs is not None:
        try:
            n_time_obs = len(time_obs)
        except Exception:
            n_time_obs = 1
        iris_response["AREA_SG"] = np.zeros(iris_response["AREA_SG"].shape)
        iris_response["AREA_SJI"] = np.zeros(iris_response["AREA_SJI"].shape)

        # 1. FUV SG effective areas
        lambran_fuv = np.array([[133.1, 135.9],[138.8, 140.8]])
        # Rough SG spectral ranges.  Setting effective area to 0 outside of these.
        shp_fuv = iris_response["COEFFS_FUV"].shape
        # Time-dependent response for shp_0[0] = 3 wavelengths
        iris_fit_fuv = np.zeros((n_time_obs, shp_fuv[0]))
        detector_fuv = "FUV"
        for j in range(shp_fuv[0]):
            iris_fit_fuv[:, j] = fit_iris_xput(time_obs, iris_response["C_F_TIME"], iris_response["COEFFS_FUV"][j, :, :])
        # Interpolate onto lambda grid, separately for each of the two FUV CCD's.
        for j in range(2):
            w_fuv = np.logical_and(iris_response["LAMBDA"].value >= lambran_fuv[j, 0], iris_response["LAMBDA"].value  <= lambran_fuv[j, 1])
            for k in range(n_time_obs):
                interpol_fuv = scipy.interpolate.interp1d(iris_response["C_F_LAMBDA"][j:j+2], np.squeeze(iris_fit_fuv[k, j:j+2]), fill_value='extrapolate')
                iris_response["AREA_SG"][0, w_fuv] = interpol_fuv(iris_response["LAMBDA"][w_fuv])

        # 2. NUV SG effective areas
        lambran_nuv = np.array([278.2, 283.5])
        # Rough SG spectral ranges.  Setting effective area to 0 outside of these.
        shp_nuv = iris_response["COEFFS_NUV"].shape
        # Time-dependent response for shp_1[0] wavelengths
        iris_fit_nuv = np.zeros((n_time_obs, shp_nuv[0]))
        detector_nuv = "NUV"
        for j in range(shp_nuv[0]):
            iris_fit_nuv[:, j] = fit_iris_xput(time_obs, iris_response["C_N_TIME"], iris_response["COEFFS_NUV"][j, :, :])
        # Interpolate onto lambda grid
        w_nuv = np.where(np.logical_and(iris_response["LAMBDA"].value >= lambran_nuv[0], iris_response["LAMBDA"].value <= lambran_nuv[1]))
        if int(iris_response["VERSION"]) <= 3:
            for k in range(n_time_obs):
                interpol_nuv =  scipy.interpolate.interp1d(iris_response["C_N_LAMBDA"][:], np.squeeze(iris_fit_nuv[k, :]), fill_value='extrapolate')
                iris_response["AREA_SG"][1, w_nuv] = interpol_nuv(iris_response["LAMBDA"][w_nuv])
        else:
            for k in range(n_time_obs):
                interpol_nuv = scipy.interpolate.CubicSpline(iris_response["C_N_LAMBDA"][:], np.squeeze(iris_fit_nuv[k, :]),  extrapolate=True, bc_type="natural", axis=0)
                iris_response["AREA_SG"][1, w_nuv] = interpol_nuv(iris_response["LAMBDA"][w_nuv])

        # 3. SJI effective areas
        if int(iris_response["VERSION"]) <= 3:  # Meaning for version 3 only
            shp_sji = iris_response["COEFFS_SJI"].shape
            for j in range(shp_sji[0]):
                # Calculate pre-launch area from the individual elements
                prelaunch_area = iris_response["GEOM_AREA"]
                for k in range(len(iris_response["INDEX_EL_SJI"][j, :])):
                    index_values0 = iris_response["INDEX_EL_SJI"][j, k]
                    prelaunch_area = prelaunch_area * iris_response["ELEMENTS"][index_values0].trans
                # Time dependent response
                iris_fit_sji = fit_iris_xput(time_obs, iris_response["C_S_TIME"][j, :, :], iris_response["COEFFS_SJI"][j, :, :])
                # Time dependetn profiles
                for k in range(n_time_obs):
                    iris_response["AREA_SJI"][j, :] = prelaunch_area * iris_fit_sji[k]
        else:  # For version 4 and above
            for nuv in range(2):
                # Calculate baseline SJI area curves
                area_sji = iris_response["GEOM_AREA"]
                for m in range(len(iris_response["INDEX_EL_SJI"][nuv*2, :])):
                    index_values1 = iris_response["INDEX_EL_SJI"][nuv*2: nuv*2+2, m]
                    area_sji = area_sji * iris_response["ELEMENTS"][index_values1].trans
                # Apply time dependent profile shape adjustment to FUV SJI
                if nuv == 0:
                    # FUV: apply FUV SG "slant", then normalize so that a weighted (2.4:1) sum at C II and Si IV gives constant response
                    weight = np.array([2.4, 1.0])  # Typical solar ratio CII : SiIV
                    wavelength = iris_response["C_F_LAMBDA"]
                    n_wavelength = len(wavelength)
                    wavelength = np.array([wavelength[0].value, (wavelength[n_wavelength-2].value*2.0 + wavelength[n_wavelength-1].value)/3.0])  # 2 wvlngts in nm
                    # Calculate baseline SG area for scaling purposes
                    area_sg = iris_response["GEOM_AREA"]
                    for n in range(len(iris_response["INDEX_EL_SG"][nuv, :])):
                        index_values2 = iris_response["INDEX_EL_SG"][nuv, n]
                        area_sg = area_sg * iris_response["ELEMENTS"][index_values2].trans
                    # SG and SJI areas at wavelength
                    interpol_sg = scipy.interpolate.interp1d(iris_response["LAMBDA"], np.squeeze(area_sg), fill_value='extrapolate')
                    area_sg2 = interpol_sg(wavelength)
                    area_sj2 = np.zeros((2, 2))
                    for n in range(2):
                        interpol_sji = scipy.interpolate.interp1d(iris_response["LAMBDA"], np.squeeze(area_sji[n]), fill_value='extrapolate')
                        area_sj2[n, :] = interpol_sji(wavelength)
                    # Calculate the normalized slant function scal, apply to asji
                    for n in range(n_time_obs):
                        # Best-estimate slant, i.e., eff.area @ wavelength / baseline SG @ wavelength
                        interpol_sg2 = scipy.interpolate.interp1d(iris_response["LAMBDA"], np.squeeze(iris_response["AREA_SG"][0, :]), fill_value='extrapolate')
                        sca2 = interpol_sg2(wavelength) / area_sg2
                        # Normalize slant so that total(wei*asj2*sca2)/total(wei*asj2)=1
                        for m in range(2):
                            sca2n = sca2 * np.sum(weight*area_sj2[m, :]) / np.sum(weight * area_sj2[m, :] * sca2)
                            interpol_sca = scipy.interpolate.interp1d(wavelength, np.squeeze(sca2n), fill_value='extrapolate')
                            sca1n = interpol_sca(iris_response["LAMBDA"])
                            sca1n = np.clip(sca1n, a_min=0, a_max=None)
                            iris_response["AREA_SJI"][m] = area_sji[m] * sca1n
                else:
                    # NUV: essentially same calculation as r.version=3
                    for n in range(n_time_obs):
                        iris_response["AREA_SJI"] = [Quantity(x, unit=u.cm**2) for x in iris_response["AREA_SJI"]]
                        area_sji = [x for x in area_sji]
                        iris_response["AREA_SJI"][2:4] = area_sji[:]
            for j in range(4):
                # SJI specific time dependency
                iris_fit_sji = fit_iris_xput(time_obs, iris_response["C_S_TIME"][j, :, :], iris_response["COEFFS_SJI"][j, :, :])
                for k in range(n_time_obs):
                    iris_response["AREA_SJI"][j] = iris_response["AREA_SJI"][j] * iris_fit_sji[k]

    if not isinstance(iris_response["AREA_SG"], Quantity):
        iris_response["AREA_SG"] = Quantity(iris_response["AREA_SG"], unit=u.cm**2)
    if not isinstance(iris_response["AREA_SJI"], Quantity):
        iris_response["AREA_SJI"] = Quantity(iris_response["AREA_SJI"], unit=u.cm**2)

    return iris_response


def fit_iris_xput(time_obs, time_cal_coeffs, cal_coeffs):
    """
    To calculate the coefficients of best-fit time function for throughput,
    for which there are two modes:
    1. Perform fit: supply xput and single element ``cal_coeffs``.
    2. Apply fit: supply full ``cal_coeffs``.

    The procedure involved in this function is as follows:
    1. The time difference (in years) is computed from the time_obs and time_cal_coeffs.
    2. A least-sqaures fit is performed to determine the best fit for the time-dependent
    effective areas given the time difference.

    Parameters
    ----------
    time_obs: a `numpy.array`
        A set of observation times as `astropy.time.Time` objects contained
        in a numpy array.

    time_cal_coeffs: a numpy array of floats (with exactly two columns)
    - Start and end times of intervals of constant ``cal_coeffs[i]``.

    cal_coeffs: a numpy array of floats (with at least two columns)
    - Coefficients of best-fit function.

    Returns
    -------
    fit_out: `numpy.array`
        Yields the fit used to compute the effective area using the input
    times ``time_obs``.

    """
    time_obs = np.array([parse_time(time_obs, format='utime')])

    size_time_cal_coeffs = time_cal_coeffs.shape
    size_cal_coeffs = cal_coeffs.shape

    if size_time_cal_coeffs[1] != 2 or size_cal_coeffs[1] < 2:
        # Raise ValueError as time coefficient have the wrong format.
        raise ValueError("Incorrect number of elements either in time_cal_coeffs or in cal_coeffs.")

    # Some time transformations.
    # Convert the time_cal_coeffs given in the .geny file into a ``astropy.time.Time``
    # object called t_cal_coeffs, so that the time differences will be in days...
    t_cal_coeffs_flat = time_cal_coeffs.flatten()
    t_cal_coeffs = [parse_time(x, format = 'utime') for x in t_cal_coeffs_flat]
    t_cal_coeffs = np.array(t_cal_coeffs).reshape(size_time_cal_coeffs)

    # Exponent for transition between exp.decay intervals.
    transition_exp = 1.5

    # For loop for carrying out the least-squares fit and computation of fit output.
    fit_out = np.zeros(len(time_obs))

    for i, t in enumerate(time_obs):
        aux_cal_coeffs = np.zeros(2*size_time_cal_coeffs[0])

        fit_out = np.zeros(len(list(time_obs)), dtype=np.float64)

        # Looking for the closest time in the calibration time intervals.
        # Differences are given in years before passing to the next stage.
        t_diff = t - t_cal_coeffs

        t_diff = t_diff.flatten()

        # To Convert to an array, quantities need to be dimensionless, hence dividing out the unit.
        t_diff = np.array([x.to(u.year).value for x in t_diff])
        w = np.where(t_diff < 0)[0][0]

        # If the t_obs is betweeen the calibration time intervals of a
        # calibration file (w % !=0) then the aux_coeffs are given by an
        # exponential (coefficient and exponential value).
        # If the t_obs is between the end calibration time interval of
        # a calibration file (cal_file_t) and the beginning calibration
        # time interval of the next calibration file (cal_file_t+1)
        # (w % 2 == 0) then, the aux_coeffs are given by 4 values
        # corresponding to a partial exponential obtained from
        # cal_file_t and a complemetary exponential obtained from the
        # cal_file_t+1
        if w % 2 != 0:  # I.e., if w is not even...
            dtt_0 = 1.
            exp_0 = np.exp(cal_coeffs[w//2, 2]*(t_diff[w-1]))
            aux_cal_coeffs[w-1:w+1] = np.array([dtt_0, dtt_0*exp_0])
        else:
            dtt_1 = (t_diff[w-1] / (t_diff[w-1] - t_diff[w]))**transition_exp
            dtt_0 = 1. - dtt_1
            exp_0 = np.exp(cal_coeffs[(w//2)-1, 2]*(t_diff[w-2]))
            exp_1 = np.exp(cal_coeffs[w//2, 2]*(t_diff[w]))
            aux_cal_coeffs[w-2:w+2] = np.array([dtt_0, dtt_0*exp_0, dtt_1, dtt_1*exp_1])
        # print(aux_cal_coeffs)
        # print(cal_coeffs[:,:2].reshape((aux_cal_coeffs.shape[0])))
        fit_out[i] = np.matmul(aux_cal_coeffs, cal_coeffs[:, :2].reshape((aux_cal_coeffs.shape[0])))

    return fit_out


@custom_model
def _gaussian1d_on_linear_bg(x, amplitude=None, mean=None, standard_deviation=None,
                             constant_term=None, linear_term=None):
    return amplitude * np.exp(-((x - mean) / standard_deviation) ** 2) + constant_term + linear_term * x


def get_detector_type(meta):
    """
    Gets the IRIS detector type from a meta dictionary.

    In this function, FUV1 and FUV2 are just assigned as FUV.

    Parameters
    ----------
    meta: dict-like
        Dictionary-like object containing entry for "detector type"

    Returns
    -------
    detector_type: `str`
       Detector type.

    """
    if "FUV" in meta["detector type"]:
        detector_type = "FUV"
    else:
        detector_type = meta["detector type"]
    return detector_type


def convert_between_DN_and_photons(old_data_arrays, old_unit, new_unit):
    """Converts arrays from IRIS DN to photons or vice versa.

    In this function, an inverse time component due to exposure time
    correction is ignored during calculations but preserved in final unit.

    Parameters
    ----------
    old_data_arrays: iterable of `numpy.ndarray`s
        Arrays of data to be converted.

    old_unit: `astropy.unit.Unit`
        Unit of data arrays.

    new_unit: `astropy.unit.Unit`
        Unit to convert data arrays to.

    Returns
    -------
    new_data_arrays: `list` of `numpy.ndarray`s
        Data arrays converted to new_unit.

    new_unit_time_accounted: `astropy.unit.Unit`
        Unit of new data arrays with any inverse time component preserved.

    """
    if old_unit == new_unit or old_unit == new_unit / u.s:
        new_data_arrays = [data for data in old_data_arrays]
        new_unit_time_accounted = old_unit
    else:
        # During calculations, the time component due to exposure
        # time correction, if it has been applied, is ignored.
        # Check here whether the time correction is present in the
        # original unit so that is carried through to new unit.
        if u.s not in (old_unit * u.s).decompose().bases:
            old_unit_without_time = old_unit * u.s
            new_unit_time_accounted = new_unit / u.s
        else:
            old_unit_without_time = old_unit
            new_unit_time_accounted = new_unit
        # Convert data and uncertainty to new unit.
        new_data_arrays = [(data * old_unit_without_time).to(new_unit).value
                           for data in old_data_arrays]
    return new_data_arrays, new_unit_time_accounted


def convert_or_undo_photons_per_sec_to_radiance(data_quantities,
        time_obs, response_version, obs_wavelength, detector_type,
        spectral_dispersion_per_pixel, solid_angle, undo=False):
    """
    Converts data quantities from counts/s to radiance (or vice versa).

    Parameters
    ----------
    data_quantities: iterable of `astropy.units.Quantity`s
        Quantities to be converted.  Must have units of counts/s or
        radiance equivalent counts, e.g. erg / cm**2 / s / sr / Angstrom.

    time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2
        Observation times of the datapoints.
        Must be in the format of, e.g.,
        time_obs parse_time('2013-09-03', format='utime'),
        which yields 1094169600.0 seconds in value.
        The argument time_obs is ignored for versions 1 and 2.

    response_version : `int`
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=4.

    obs_wavelength: `astropy.units.Quantity`
        Wavelength at each element along spectral axis of data quantities.

    detector_type: `str`
        Detector type: 'FUV', 'NUV', or 'SJI'.

    spectral_dispersion_per_pixel: scalar `astropy.units.Quantity`
        spectral dispersion (wavelength width) of a pixel.

    solid_angle: scalar `astropy.units.Quantity`
        Solid angle corresponding to a pixel.

    undo: `bool`
        If False, converts counts/s to radiance.
        If True, converts radiance to counts/s.
        Default=False

    Returns
    -------
    new_data_quantities: `list` of `astropy.units.Quantity`s
        Data quantities converted to radiance or counts/s
        depending on value of undo kwarg.

    """
    # Check data quantities are in the right units.
    if undo is True:
        for i, data in enumerate(data_quantities):
            if not data.unit.is_equivalent(RADIANCE_UNIT):
                raise ValueError(
                    "Invalid unit provided.  As kwarg undo=True, "
                    "unit must be equivalent to {0}.  Error found for {1}th element "
                    "of data_quantities. Unit: {2}".format(RADIANCE_UNIT, i, data.unit))
    else:
        for data in data_quantities:
            if data.unit != u.photon/u.s:
                raise ValueError(
                    "Invalid unit provided.  As kwarg undo=False, "
                    "unit must be equivalent to {0}.  Error found for {1}th element "
                    "of data_quantities. Unit: {2}".format(u.photon/u.s, i, data.unit))
    photons_per_sec_to_radiance_factor = calculate_photons_per_sec_to_radiance_factor(time_obs,
        response_version, obs_wavelength, detector_type, spectral_dispersion_per_pixel, solid_angle)
    # Change shape of arrays so they are compatible for broadcasting
    # with data and uncertainty arrays.
    photons_per_sec_to_radiance_factor = \
    _reshape_1D_wavelength_dimensions_for_broadcast(photons_per_sec_to_radiance_factor, data_quantities[0].ndim)
    # Perform (or undo) radiometric conversion.
    if undo is True:
        new_data_quantities = [(data / photons_per_sec_to_radiance_factor).to(u.photon/u.s)
                               for data in data_quantities]
    else:
        new_data_quantities = [(data*photons_per_sec_to_radiance_factor).to(RADIANCE_UNIT)
                               for data in data_quantities]
    return new_data_quantities


def calculate_photons_per_sec_to_radiance_factor(time_obs, response_version,
        wavelength, detector_type, spectral_dispersion_per_pixel, solid_angle):
    """
    Calculates multiplicative factor that converts counts/s to radiance for given wavelengths.

    Parameters
    ----------
    time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2
        Observation times of the datapoints.
        Must be in the format of, e.g.,
        time_obs=parse_time('2013-09-03', format='utime'),
        which yields 1094169600.0 seconds in value.
        The argument time_obs is ignored for versions 1 and 2.

    response_version : `int`
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=4.

    wavelength: `astropy.units.Quantity`
        Wavelengths for which counts/s-to-radiance factor is to be calculated

    detector_type: `str`
        Detector type: 'FUV' or 'NUV'.

    spectral_dispersion_per_pixel: scalar `astropy.units.Quantity`
        spectral dispersion (wavelength width) of a pixel.

    solid_angle: scalar `astropy.units.Quantity`
        Solid angle corresponding to a pixel.

    Returns
    -------
    radiance_factor: `astropy.units.Quantity`
        Mutliplicative conversion factor from counts/s to radiance units
        for input wavelengths.

    """
    # Get effective area and interpolate to observed wavelength grid.
    eff_area_interp = _get_interpolated_effective_area(time_obs,
        response_version, detector_type, obs_wavelength=wavelength)
    # Return radiometric conversed data assuming input data is in units of photons/s.
    return constants.h * constants.c / wavelength / u.photon / spectral_dispersion_per_pixel / eff_area_interp / solid_angle


def _get_interpolated_effective_area(time_obs, response_version, detector_type, obs_wavelength, *args, **kwargs):
    """
    To compute the interpolated time-dependent effective area.

    Parameters
    ----------
    time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2
        Observation times of the datapoints.
        Must be in the format of, e.g.,
        time_obs parse_time('2013-09-03', format='utime'),
        which yields 1094169600.0 seconds in value.
        The argument time_obs is ignored for versions 1 and 2.

    response_version : `int`
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=4.

    detector_type: `str`
        Detector type: 'FUV' or 'NUV'.

    obs_wavelength: `astropy.units.Quantity`
        The wavelength at which the observation has been taken in Angstroms.

    Returns
    -------
    eff_area_interp: `numpy.array`
        The effective area(s) determined by interpolation with a spline fit.

    """
    # Generalizing to the time of obs.
    time_obs = time_obs
    response_version = response_version
    iris_response = get_iris_response(time_obs, response_version,  *args, **kwargs)
    if detector_type == "FUV":
        detector_type_index = 0
    elif detector_type == "NUV":
        detector_type_index = 1
    else:
        raise ValueError("Detector type not recognized.")
    eff_area = iris_response["AREA_SG"][detector_type_index, :]
    response_wavelength = iris_response["LAMBDA"]
    # Interpolate the effective areas to cover the wavelengths
    # at which the data is recorded:
    eff_area_interp_base_unit = u.Angstrom
    tck = interpolate.splrep(response_wavelength.to(eff_area_interp_base_unit).value,
                             eff_area.to(eff_area_interp_base_unit ** 2).value, s=0)
    eff_area_interp = interpolate.splev(
        obs_wavelength.to(eff_area_interp_base_unit).value, tck) * (eff_area_interp_base_unit ** 2)
    return eff_area_interp


def _reshape_1D_wavelength_dimensions_for_broadcast(wavelength, n_data_dim):
    if n_data_dim == 1:
        pass
    elif n_data_dim == 2:
        wavelength = wavelength[np.newaxis, :]
    elif n_data_dim == 3:
        wavelength = wavelength[np.newaxis, np.newaxis, :]
    else:
        raise ValueError("IRISSpectrogram dimensions must be 2 or 3.")
    return wavelength
