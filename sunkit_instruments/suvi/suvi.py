import os
from os.path import expanduser
import matplotlib.pyplot as plt
import numpy
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy import interpolate
from scipy.ndimage import gaussian_filter
import sunpy.map
from astropy import units as u
from astropy.io import fits
from astropy.io.fits.verify import VerifyError
from astropy.time import Time, TimeDelta
from sunpy.util.exceptions import warn_user
import sunkit_instruments

SOLAR_CLASSES = [('unlabeled', 0),
                 ('outer_space', 1),
                 ('bright_region', 3),
                 ('filament', 4),
                 ('prominence', 5),
                 ('coronal_hole', 6),
                 ('quiet_sun', 7),
                 ('limb', 8),
                 ('flare', 9)]

SOLAR_CLASS_NAME = {number: theme for theme, number in SOLAR_CLASSES}

SOLAR_COLORS = {"unlabeled": "white",
                "outer_space": "black",
                "bright_region": "#F0E442",
                "filament": "#D55E00",
                "prominence": "#E69F00",
                "coronal_hole": "#009E73",
                "quiet_sun": "#0072B2",
                "limb": "#56B4E9",
                "flare": "#CC79A7"}


__all__ = ["despike_L1b_image",
           "files_to_map",
           "get_response",
           "plot_thematic_map"]


def despike_L1b_image(the_input, filter_width=7, return_map=False):
    """
    Despike SUVI L1b data and return either a despiked
    `~numpy.ndarray` or a despiked `~sunpy.map.Map`.

    .. note::
        If the the_input is a string: the type of file is determined
        by pattern matching in the filenames, e.g. "-L1b-Fe171"
        for a 171 L1b file. If this pattern is not found in the
        filename, the file will not be recognized.

    .. note::
        The despiking relies on the presence of the data quality
        flags (DQF) in the first extension of a SUVI L1b FITS file.
        Early in the mission, the DQF extension was not present
        yet, so the despiking cannot be done with this function
        for those early files.

    Parameters
    ----------
    the_input: `str` or `tuple` of 2-D `~numpy.ndarray` containing (data, DQF)
        Filename or tuple of (data, DQF) to despike.

    filter_width: `int`, optional. Default: 7.
        The filter width for the Gaussian filter. If NaNs are still
        present in the despiked image, try increasing this value.

    return_map: `bool`, optional. Default: False.
        If True, returns a `~sunpy.map.Map` instead of a `~numpy.ndarray`.
        Ignored if the_input is not a filename string.

    Returns
    -------
    despiked_image: `~numpy.ndarray`
        The despiked L1b image as a numpy array.

    map: `~sunpy.map.Map` if return_map=True.
        The despiked L1b image as a `~sunpy.map.Map`.
    """

    if isinstance(the_input, tuple):
        image = the_input[0]
        dqf_mask = the_input[1]
    elif isinstance(the_input, str):
        L1b_matches = ['-L1b-Fe093', '-L1b-Fe131', '-L1b-Fe171', '-L1b-Fe195', '-L1b-Fe284', '-L1b-He303']
        if any(fn in os.path.basename(the_input) for fn in L1b_matches):
            header, image, dqf_mask = sunkit_instruments.suvi.read_suvi(the_input, return_DQF=True)
        else:
            raise ValueError("File "+the_input+" does not look like a SUVI L1b file.")
    else:
        raise TypeError("Input must be string or tuple of 2-D arrays.")

    image_with_nans = numpy.copy(image)
    image_with_nans[numpy.where(dqf_mask == 4)] = numpy.nan
    indices = numpy.where(numpy.isnan(image_with_nans))
    image_gaussian_filtered = gaussian_filter(image, filter_width)
    despiked_image = numpy.copy(image_with_nans)
    despiked_image[indices] = image_gaussian_filtered[indices]

    if return_map and isinstance(the_input, str):
        return sunpy.map.Map(despiked_image, header)
    else:
        return despiked_image


def files_to_map(files, sort_files=True, verbose=False, despike_L1b=False,
                 only_long_exposures=False, only_short_exposures=False, only_short_flare_exposures=False):
    """
    Read SUVI L1b FITS or netCDF files or L2 HDR composite FITS files and
    return a `~sunpy.map.Map` or a `~sunpy.map.MapSequence`. For SUVI L1b
    FITS files, the broken FITS header is fixed automatically (broken
    because of the wrong implementation of the CONTINUE convention).

    .. note::
        The first file in the (sorted, if sort_files=True) list determines what
        will be accepted further on, i.e. L2 HDR composites or L1b files. If L1b
        files are appearing in a file list that started with an L2 HDR composite,
        they will be rejected (and vice versa). The type of file is determined
        by pattern matching in the filenames, e.g. "-L1b-Fe171" for a 171 L1b file
        and "-l2-ci171" for a 171 L2 HDR composite. If those patterns are not found
        in the filename, the files will not be recognized.

    Parameters
    ----------
    files: `str` or `list` of `str`
        File(s) to read.

    sort_files: `bool`, optional. Default: True.
        If True, sorts the input file list (ascending).

    verbose: `bool`, optional. Default: False.
        If True, prints the filenames while reading.

    despike_L1b: `bool`, optional. Default: False.
        If True and input is L1b, data will get despiked
        with the standard filter_width=7. Can not be used
        for early SUVI files where the DQF extension is
        missing.

    only_long_exposures: `bool`, optional. Default: False.
        If True, only long exposure L1b files from the input list will be
        accepted and converted to a map. Ignored for L2 HDR composites.

    only_short_exposures: `bool`, optional. Default: False.
        If True, only short exposure L1b files from the input list will be
        accepted and converted to a map. Ignored for L2 HDR composites and
        any wavelengths other than 94 and 131 (because for everything >131,
        there are no observations that are labeled "short", only "long" and
        "short_flare").

    only_short_flare_exposures: `bool`, optional. Default: False.
        If True, only short flare exposure L1b files from the input list will
        be accepted and converted to a map. Ignored for L2 HDR composites.

    Returns
    -------
    `~sunpy.map.Map`, `~sunpy.map.MapSequence`, or `None`.
        A map (sequence) of the SUVI data, or None if no
        data was found matching the given criteria.
    """

    # If it is just one filename as a string, convert it to a list.
    if isinstance(files, str):
        files = [files]
    else:
        if sort_files:
            files = sorted(files)

    # Based on the filename of the first file, determine if we are
    # dealing with L1b files or HDR composites (or neither).
    composite_matches = ['-l2-ci094', '-l2-ci131', '-l2-ci171', '-l2-ci195', '-l2-ci284', '-l2-ci304']
    L1b_matches = ['-L1b-Fe093', '-L1b-Fe131', '-L1b-Fe171', '-L1b-Fe195', '-L1b-Fe284', '-L1b-He303']
    if any(fn in os.path.basename(files[0]) for fn in composite_matches):
        composites = True
    elif any(fn in os.path.basename(files[0]) for fn in L1b_matches):
        composites = False
    else:
        raise ValueError("First file "+files[0]+" does not look like a SUVI L1b file or L2 HDR composite.")

    datas = []
    headers = []
    for tmp_file in files:
        if verbose:
            print('Reading', tmp_file)
        # Test for L1b or composite based on the filename
        if composites:
            if any(fn in os.path.basename(tmp_file) for fn in composite_matches):
                tmp_header, tmp_data, _ = sunkit_instruments.suvi.read_suvi(tmp_file)
                datas.append(tmp_data)
                headers.append(tmp_header)
            else:
                warn_user("File "+tmp_file+" does not look like a SUVI L2 HDR composite. Skipping.")
        else:
            if any(fn in os.path.basename(tmp_file) for fn in L1b_matches):
                if despike_L1b:
                    tmp_header, tmp_data, dqf_mask = sunkit_instruments.suvi.read_suvi(tmp_file, return_DQF=True)
                    tmp_data = despike_L1b_image((tmp_data, dqf_mask))
                else:
                    tmp_header, tmp_data, _ = sunkit_instruments.suvi.read_suvi(tmp_file)
                if only_long_exposures:
                    if 'long_exposure' in tmp_header['SCI_OBJ']:
                        datas.append(tmp_data)
                        headers.append(tmp_header)
                elif only_short_exposures:
                    if 'short_exposure' in tmp_header['SCI_OBJ']:
                        datas.append(tmp_data)
                        headers.append(tmp_header)
                elif only_short_flare_exposures:
                    if 'short_flare_exposure' in tmp_header['SCI_OBJ']:
                        datas.append(tmp_data)
                        headers.append(tmp_header)
                else:
                    datas.append(tmp_data)
                    headers.append(tmp_header)
            else:
                warn_user("File "+tmp_file+" does not look like a SUVI L1b file. Skipping.")

    if len(datas) == 1 and len(headers) == 1:
        # Make a single map if it is just one file
        suvimap = sunpy.map.Map(datas[0], headers[0])
        return suvimap
    elif len(datas) > 1 and len(headers) > 1:
        # Make a map sequence for multiple files
        suvimap = sunpy.map.Map(list(zip(datas, headers)), sequence=True)
        return suvimap
    else:
        warn_user("List of data/headers is empty.")
        return None


def get_response(the_input, spacecraft=16, ccd_temperature=-60., exposure_type='long'):
    """
    Get the SUVI instrument response for a specific wavelength channel,
    spacecraft, CCD temperature, and exposure type. the_input can either
    be an L1b filename (FITS or netCDF), in which case all of those
    parameters are read automatically from the metadata, or the parameters
    can be passed manually, with the_input specifying the desired wavelength
    channel.

    Parameters
    ----------
    the_input: `str` or `int`.
        Either an L1b filename (FITS or netCDF), or an integer specifying the
        wavelength channel. Valid wavelength channels: 94, 131, 171, 195, 284, and 304.

    spacecraft: `str` or `int`, optional. Default: 16.
        Which spacecraft. Can either be the full spacecraft name, or just the number.
        Valid: "GOES-16", "GOES-17", 16, and 17.

    ccd_temperature: `float`, optional. Default: -60.
        The CCD temperature, in degrees Celsius. Needed for getting
        the correct gain number.

    exposure_type: `str`, optional. Default: "long".
        The exposure type of the SUVI image. Can be:
        "long", "short", "short_flare" for 94 and 131;
        "long", "short_flare" for 171, 195, 284, and 304.
        The exposure type is needed for the correct focal plane
        filter selection.

    Returns
    -------
    `dict` with the instrument response information.
        Keys: "wavelength", "effective_area", "response",
        "wavelength_channel", "spacecraft", "ccd_temperature",
        "exposure_type", "flight_model", "gain", "solid_angle",
        "geometric_area", "filter_setup".
    """

    if isinstance(the_input, str):
        hdr, _, _ = sunkit_instruments.suvi.read_suvi(the_input, return_header_only=True)
        wavelength_channel = int(hdr['WAVELNTH'])
        spacecraft = int(hdr['TELESCOP'].replace(' ','').replace('G',''))
        ccd_temperature = (hdr['CCD_TMP1'] + hdr['CCD_TMP2'])/2.
        exposure_type = '_'.join(hdr['SCI_OBJ'].replace(" ","").split(sep='_')[3:]).replace('_exposure','')
    elif isinstance(the_input, int):
        wavelength_channel = the_input
    else:
        raise TypeError("Input not recognized, must be str for filename or int for wavelength channel.")

    valid_wavelength_channels = [94, 131, 171, 195, 284, 304]
    if wavelength_channel not in valid_wavelength_channels:
        raise ValueError("Wavelength channel"+str(wavelength_channel)+" not recognized. Valid: "+str(valid_wavelength_channels))

    valid_spacecraft = ['GOES-16', 'GOES-17', 16, 17]
    if spacecraft not in valid_spacecraft:
        raise ValueError("Spacecraft "+str(spacecraft)+" not recognized. Valid: "+str(valid_spacecraft))

    if isinstance(spacecraft, str):
        spacecraft = int(spacecraft.split('-')[1])

    flight_model = {16: 'FM1', 17: 'FM2'}

    # The setup for filter wheel 1 and 2 for the different exposure types
    filter_setup = { 94: {'long':        {'FW1': 'thin_zirconium', 'FW2': 'open'},
                          'short':       {'FW1': 'thin_zirconium', 'FW2': 'open'},
                          'short_flare': {'FW1': 'thin_zirconium', 'FW2': 'thin_zirconium'}},
                    131: {'long':        {'FW1': 'thin_zirconium', 'FW2': 'open'},
                          'short':       {'FW1': 'thin_zirconium', 'FW2': 'open'},
                          'short_flare': {'FW1': 'thin_zirconium', 'FW2': 'thin_zirconium'}},
                    171: {'long':        {'FW1': 'thin_aluminum',  'FW2': 'open'},
                          'short_flare': {'FW1': 'thin_aluminum',  'FW2': 'thin_aluminum'}},
                    195: {'long':        {'FW1': 'thin_aluminum',  'FW2': 'open'},
                          'short_flare': {'FW1': 'thin_aluminum',  'FW2': 'thin_aluminum'}},
                    284: {'long':        {'FW1': 'thin_aluminum',  'FW2': 'open'},
                          'short_flare': {'FW1': 'thin_aluminum',  'FW2': 'thin_aluminum'}},
                    304: {'long':        {'FW1': 'thin_aluminum',  'FW2': 'open'},
                          'short_flare': {'FW1': 'thin_aluminum',  'FW2': 'thin_aluminum'}}}

    path_to_files = os.path.join(os.path.dirname(sunkit_instruments.__file__), 'suvi', 'data')
    eff_area_file = os.path.join(path_to_files, 'SUVI_'+flight_model[spacecraft]+'_'+str(wavelength_channel)+'A_eff_area.txt')
    gain_file = os.path.join(path_to_files, 'SUVI_'+flight_model[spacecraft]+'_gain.txt')

    # Get effective areas
    eff_area = numpy.loadtxt(eff_area_file, skiprows=12)
    wave = eff_area[:,0] * u.Angstrom
    if filter_setup[wavelength_channel][exposure_type]['FW2'] == 'open':
        effective_area = eff_area[:,1] * u.cm * u.cm
    else:
        effective_area = eff_area[:,2] * u.cm * u.cm

    # Get gain
    gain_table = numpy.loadtxt(gain_file, skiprows=7)
    temp_x = gain_table[:,0]
    gain_y = gain_table[:,1]
    gain_vs_temp = interpolate.interp1d(temp_x, gain_y)
    gain = gain_vs_temp(ccd_temperature)

    geometric_area = 19.362316 * u.cm * u.cm
    solid_angle = ((2.5/3600. * (numpy.pi/180.))**2.) * u.sr
    master_e_per_phot = ((6.626068e-34 * (u.J/u.Hz)) * (2.99792458e8 * (u.m/u.s)))/(wave.to(u.m) * ((u.eV.to(u.J, 3.65)) * u.J))

    response = effective_area * (master_e_per_phot/gain)

    response_info = {'wavelength': wave, 'effective_area': effective_area * (u.ct/u.ph), 'response': response,
                     'wavelength_channel': wavelength_channel, 'spacecraft': 'GOES-'+str(spacecraft),
                     'ccd_temperature': ccd_temperature * u.deg_C, 'exposure_type': exposure_type,
                     'flight_model': flight_model[spacecraft], 'gain': numpy.float64(gain),
                     'solid_angle': solid_angle, 'geometric_area': geometric_area,
                     'filter_setup': filter_setup[wavelength_channel][exposure_type]}

    return response_info


def plot_thematic_map(input_filename, timestamp=True, legend=True, figsize=(10,10)):
    """
    Read a SUVI L2 Thematic Map FITS file and plot it.

    .. note::
        SUVI L2 Thematic Maps are recognized by pattern matching in the
        filenames, i.e. it must contain "-l2-thmap". If this pattern is
        not found in the filename, it will not be recognized.

    .. note::
       GOES-16 L2 Thematic Maps are available
       `here. <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/suvi-l2-thmap/>`__

    Parameters
    ----------
    input_filename: `str`
        File to read.

    timestamp: `bool`, optional. Default: True.
        If True, plots the timestamp of the observation.

    legend: `bool`, optional. Default: True.
        If True, plots a legend with the different classes.

    figsize: `tuple`, optional. Default: (10, 10).
        The size of the figure. Should not be smaller than
        (7, 7) to look decent.
    """

    if '-l2-thmap' in os.path.basename(input_filename):
        hdu = fits.open(input_filename)
        thmap_data = hdu[0].data
        if timestamp:
            time_st = hdu[0].header['DATE-OBS'][0:19]
        hdu.close()

        colortable = [SOLAR_COLORS[SOLAR_CLASS_NAME[i]] if i in SOLAR_CLASS_NAME else 'black'
                      for i in range(max(list(SOLAR_CLASS_NAME.keys())) + 1)]
        cmap = ListedColormap(colortable)

        fig, ax = plt.subplots(constrained_layout=True, figsize=figsize)
        ax.imshow(thmap_data, origin='lower', cmap=cmap, vmin=-1, vmax=len(colortable), interpolation='none')
        ax.set_axis_off()
        if timestamp:
            ax.text(10, 1245, time_st, fontsize=14, color='white')
        if legend:
            legend_elements = [Patch(facecolor=color, edgecolor="black", label=label.replace("_", " "))
                               for label, color in SOLAR_COLORS.items()]
            ax.legend(handles=legend_elements, loc='upper right', ncol=3, fancybox=True, shadow=False)
        fig.show()
    else:
        raise ValueError("File "+input_filename+" does not look like a SUVI L2 Thematic Map.")
