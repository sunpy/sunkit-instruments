import os
import gzip
import tempfile
from os.path import expanduser
import matplotlib.pyplot as plt
import numpy
import requests
from bs4 import BeautifulSoup
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
           "download_data_from_NOAA",
           "files_to_map",
           "fix_L1b_header",
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

    
# helper function: parser for the SUVI websites using BeautifulSoup
def _list_url_directory(url, ext=''):
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [url + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]


def download_data_from_NOAA(date_time, spacecraft=16, wavelength=171, outdir=None,
                            composites=False, query_only=False, verbose=False):
    """
    Download SUVI L1b or L2 HDR composite FITS files from the NOAA websites directly.
    Intended as an alternative to getting it from the VSO with fido. For L1b
    files, it will only download long exposures, no short/short flare exposures.

    .. note::
        Existing files with the same filename will be overwritten.

    .. note::
        Using the JSOC-style format below does not mean that the cadence will
        always be exactly what was requested. The code only tries to download
        the nearest file to any datetime given, i.e. if there is data missing
        for example, the code might try and download the same file over and
        over (because one file might always be the closest one). There will
        be a warning printed if that is the case.

    .. note::
        The standard observing sequence for SUVI is 4 minutes. Using 4 minutes
        (and multiples thereof) should therefore be fairly safe to use for all
        channels. 94 is observed with a higher, regular 2 min cadence.
        195 is observed with an even higher, but irregular cadence (70 s, 60 s,
        60 s, 50 s, and repeating). In other words: if higher cadence data is
        needed, it can be obtained in some channels. However, only for L1b data,
        because all L1b files get averaged into one L2 HDR composite in a 4 min
        window, and for 195 data, it can be tricky to get the timing right with
        the JSOC-style query because of the irregular cadence.

    Parameters
    ----------
    date_time: `str` or `list`.
        date_time can have the following formats:

        -A JSOC-style string with start time, timespan, and cadence:
        ``'2021-02-03T12:30:00/3h@20m'``
        Cadence is optional (defaults to 4 min if not given).
        Valid values for the time units are m, h, and d.

        -A single datetime string in isot format (without milliseconds):
        ``'2021-02-03T12:30:00'``

        -A list of several datetime strings in isot format (without milliseconds):
        ``['2021-02-03T12:30:00', '2021-04-23T11:43:00', '2021-05-11T17:05:00']``

    spacecraft: `int`, optional. Default: 16. Valid numbers: [16, 17].
        The spacecraft to get the data from.

    wavelength: `int`, optional. Default: 171. Valid numbers: [94, 131, 171, 195, 284, 304].
        The wavelength channel to get the data from.

    outdir: `str`, optional. Default: user's home directory.
        The directory to write the FITS files to.

    composites: `bool`, optional. Default: False.
        If True, the function will not download long exposure L1b files,
        but L2 HDR composite images instead.

    query_only: `bool`, optional. Default: False.
        If True, the function does not download anything, only queries the NOAA website
        and prints the filenames that would be downloaded.

    verbose: `bool`, optional. Default: False.
        If True, the user gets slightly more information about what is happening.
    """
    spacecraft_numbers = [16, 17]

    if composites:
        wvln_path = { 94:'suvi-l2-ci094', 131:'suvi-l2-ci131', 171:'suvi-l2-ci171', \
                     195:'suvi-l2-ci195', 284:'suvi-l2-ci284', 304:'suvi-l2-ci304'}
    else:
        wvln_path = { 94:'suvi-l1b-fe094', 131:'suvi-l1b-fe131', 171:'suvi-l1b-fe171', \
                     195:'suvi-l1b-fe195', 284:'suvi-l1b-fe284', 304:'suvi-l1b-he304'}

    if outdir is not None:
        if not query_only:
            if not os.path.exists(outdir):
                raise FileNotFoundError("Output directory "+outdir+" not found.")
    else:
        outdir = expanduser("~")

    # Add the path separator at the end if it is not there yet
    if outdir[-1] != os.path.sep:
        outdir = outdir + os.path.sep

    # this should stay the same for now
    baseurl1 = 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes'
    if composites:
        baseurl2 = '/l2/data/'
        ext = '.fits'
    else:
        baseurl2 = '/l1b/'
        ext = '.fits.gz'

    if spacecraft not in spacecraft_numbers:
        raise ValueError("Invalid spacecraft number: "+str(spacecraft)+". Valid values are: 16, 17.")

    if wavelength not in wvln_path:
        raise ValueError("Invalid wavelength: "+str(wavelength)+". Valid values are: 94, 131, 171, 195, 284, 304.")

    # Figure out what kind of date_time was given.
    if isinstance(date_time, str):
        # Check if it is a JSOC-style query
        if len(date_time.split('/')) == 2:
            if len(date_time.split('@')) == 2:
                cadence_string = date_time.split('@')[1]
                timespan_string = date_time.split('@')[0].split('/')[1]
                cadence = float(cadence_string[:-1])
                cadence_unit = cadence_string[-1]
                if cadence_unit == 'm':
                    cadence = cadence*60.
                elif cadence_unit == 'h':
                    cadence = cadence*60.*60.
                elif cadence_unit == 'd':
                    cadence = cadence*60.*60*24.
                else:
                    raise ValueError('Not a valid time unit (must be m, h, or d).')
            else:
                cadence = 240.
                timespan_string = date_time.split('/')[1]

            timespan = float(timespan_string[:-1])
            timespan_unit = timespan_string[-1]
            if timespan_unit == 'm':
                timespan = timespan*60.
            elif timespan_unit == 'h':
                timespan = timespan*60.*60.
            elif timespan_unit == 'd':
                timespan = timespan*60.*60*24.
            else:
                raise ValueError('Not a valid time unit (must be m, h, or d).')

            t0 = Time(date_time.split('/')[0], scale='utc', format='isot')
            tmp_timestamp = []
            counter = 0
            while counter*cadence <= timespan:
                tmp_timestamp.append(counter*cadence)
                counter += 1

            timestamp = t0+TimeDelta(tmp_timestamp, format='sec')
            urls = []
            for time in timestamp:
                urls.append(baseurl1+str(spacecraft)+baseurl2+wvln_path[wavelength]+'/'+time.value[0:10].replace('-','/')+'/')

        else:
            # Only one date, and no JSOC-style query
            timestamp = [Time(date_time, scale='utc', format='isot')]
            urls = [baseurl1+str(spacecraft)+baseurl2+wvln_path[wavelength]+'/'+date_time[0:10].replace('-','/')+'/']

    elif isinstance(date_time, list):
        # if the argument was a list of dates
        timestamp = []
        urls = []
        for this_date in date_time:
            timestamp.append(Time(this_date, scale='utc', format='isot'))
            urls.append(baseurl1+str(spacecraft)+baseurl2+wvln_path[wavelength]+'/'+this_date[0:10].replace('-','/')+'/')

    # Before we run, check if all of the websites are there.
    # Cook the urls down to unique values. To do that, convert
    # to a numpy array, use numpy.unique, and then convert back
    # to a list. Tried by using conversion to a set first,
    # but that doesn't keep the correct order for the dates.
    urls_arr = numpy.array(urls)
    urls_unique = numpy.unique(urls_arr).tolist()
    all_files  = []
    start_time = []
    if not composites:
        end_time = []
    for url in urls_unique:
        request = requests.get(url)
        if not request.status_code == 200:
            warn_user("Website not found: "+url)
            continue
        else:
            # If all of the websites were found, go ahead and make lists of files and dates.
            if verbose:
                print("Querying", url, "for SUVI files...")
            if composites:
                for file in _list_url_directory(url, ext):
                    all_files.append(file)
                    file_base = os.path.basename(file)
                    start_time.append(url[-11:-1].replace('/','-')+'T'+file_base[31:33]+':'+file_base[33:35]+':'+\
                                      file_base[35:37]+'.000')
            else:
                for file in _list_url_directory(url, ext):
                    all_files.append(file)
                    file_base = os.path.basename(file)
                    start_time.append(url[-11:-1].replace('/','-')+'T'+file_base[30:32]+':'+file_base[32:34]+':'+\
                                      file_base[34:36]+'.'+file_base[36]+'00')
                    end_time.append(url[-11:-1].replace('/','-')+'T'+file_base[46:48]+':'+file_base[48:50]+':'+\
                                    file_base[50:52]+'.'+file_base[52]+'00')

    # Make astropy time objects from the start and end times, compute the exposure time from that (only for L1b files).
    start_time = Time(start_time, scale='utc', format='isot')

    if composites:
        these_files = numpy.array(all_files)

        # Now go through all of the requested times and download/print the files.
        # Skip the download if the current file is the same one as the last.
        last_file = ' '
        for time in timestamp:
            delta_t = time-start_time
            which_file = numpy.abs(delta_t).argmin()
            if os.path.basename(these_files[which_file]) == last_file:
                if verbose:
                   warn_user("File is the same as the last one, skipping download.")
                continue
            else:
                last_file = os.path.basename(these_files[which_file])
            if query_only:
                print('Composite: ', these_files[which_file])
            else:
                if verbose:
                    print('Composite: ', these_files[which_file])
                try:
                    f = requests.get(these_files[which_file])
                    open(outdir+os.path.basename(these_files[which_file]), 'wb').write(f.content)
                finally:
                    f.close()

    else:
        end_time = Time(end_time, scale='utc', format='isot')
        exposure_time = end_time-start_time
        # Get the long exposures for the L1b files.
        long_exposures = numpy.where(numpy.around(exposure_time.sec) == 1)
        long_exposure_files = numpy.array(all_files)[long_exposures]

        # Now go through all of the requested times and download/print the files.
        # Skip the download if the current file is the same one as the last.
        last_file = ' '
        for time in timestamp:
            delta_t = time-start_time[long_exposures]
            which_file = numpy.abs(delta_t).argmin()
            if os.path.basename(long_exposure_files[which_file]) == last_file:
                if verbose:
                    warn_user("File is the same as the last one, skipping download.")
                continue
            else:
                last_file = os.path.basename(long_exposure_files[which_file])
            if query_only:
                print('Long exposure: ', long_exposure_files[which_file])
            else:
                if verbose:
                    print('Long exposure: ', long_exposure_files[which_file])
                try:
                    f = requests.get(long_exposure_files[which_file])
                    open(outdir+os.path.basename(long_exposure_files[which_file]), 'wb').write(f.content)
                finally:
                    f.close()


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
                tmp_header, tmp_data = sunkit_instruments.suvi.read_suvi(tmp_file)
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
                    tmp_header, tmp_data = sunkit_instruments.suvi.read_suvi(tmp_file)
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


def fix_L1b_header(input_filename):
    """
    Fix a SUVI L1b FITS file header (broken due to the wrong
    implementation of the CONTINUE keyword convention).

    .. note::
        astropy versions <=4.2.0 will do this faster, because we can
        still use the `astropy.io.fits.header.Header.to_string()` method.
        Starting with astropy version 4.2.1, the
        `astropy.io.fits.header.Header.to_string()` method will not work
        anymore due to FITS header consistency checks that cannot
        be overridden. The solution for that is not very elegant
        in the code here (reading the FITS file directly as bytes
        until we hit a UnicodeDecodeError), but it is the only one
        that works as far as I know.

    .. note::
        If the input file it gzipped, an unzipped file in the default
        tmp directory will be used and deleted afterwards.

    Parameters
    ----------
    input_filename: `str`
        Filename of the L1b file with the corrupt FITS header.

    Returns
    -------
    hdr_corr: `astropy.io.fits.header.Header`
        Corrected FITS header.
    """
    try:
        # First try it with the astropy .to_string() method, as this is the easiest.
        hdr = fits.getheader(input_filename)
        hdr_str = hdr.tostring()
    except VerifyError:
        # Read the file manually as bytes until we hit a UnicodeDecodeError, i.e.
        # until we reach the data part. Since astropy version 4.2.1, we can't use
        # the .to_string() method anymore because of FITS header consistency checks
        # that cannot be overridden, and they won't fix it unfortunately. If the
        # input file is a .gz file, we need to unpack it first to the tmp directory.
        temp_dir = tempfile.gettempdir()

        split_filename = os.path.splitext(os.path.basename(input_filename))
        if split_filename[1] == '.gz':
            is_gz_file = True
            with gzip.open(input_filename, 'r') as f_in, open(temp_dir + split_filename[0], 'wb') as f_out:
                f_out.write(f_in.read())
            file_to_open = temp_dir + split_filename[0]
        else:
            is_gz_file = False
            file_to_open = input_filename

        hdr_str = ''
        with open(file_to_open, 'rb') as in_file:
            counter = 1
            while True:
                try:
                    this_line = in_file.read(counter)
                    this_str = this_line.decode("utf-8")
                    hdr_str += this_str
                    counter += 1
                except UnicodeDecodeError:
                    break
        in_file.close()
        if is_gz_file:
            os.remove(file_to_open)

    # Make a list of strings with a length of 80
    hdr_list = [hdr_str[i:i+80] for i in range(0, len(hdr_str), 80)]

    # Remove all the empty entries
    while(" "*80 in hdr_list) :
        hdr_list.remove(" "*80)

    # Make a new string list where we put all the information together correctly
    hdr_list_new = []
    for count, item in enumerate(hdr_list):
        if count <= len(hdr_list)-2:
            if hdr_list[count][0:8] != 'CONTINUE' and hdr_list[count+1][0:8] != 'CONTINUE':
                hdr_list_new.append(hdr_list[count])
            else:
                if hdr_list[count][0:8] != 'CONTINUE' and hdr_list[count+1][0:8] == 'CONTINUE':
                    ampersand_pos = hdr_list[count].find('&')
                    if ampersand_pos != -1:
                        new_entry = hdr_list[count][0:ampersand_pos]
                    else:
                        # Raise exception here because there should be an ampersand at the end of a CONTINUE'd keyword
                        raise RuntimeError("There should be an ampersand at the end of a CONTINUE'd keyword.")

                    tmp_count = 1
                    while hdr_list[count+tmp_count][0:8] == 'CONTINUE':
                        ampersand_pos = hdr_list[count+tmp_count].find('&')
                        if ampersand_pos != -1:
                            first_sq_pos = hdr_list[count+tmp_count].find("'")
                            if first_sq_pos != -1:
                                new_entry = new_entry+hdr_list[count+tmp_count][first_sq_pos+1:ampersand_pos]
                            else:
                                # Raise exception here because there should be a single quote after CONTINUE
                                raise RuntimeError("There should be two single quotes after CONTINUE. Did not find any.")

                        else:
                            # If there is no ampersand at the end anymore, it means the entry ends here.
                            # Read from the first to the second single quote in this case.
                            first_sq_pos = hdr_list[count+tmp_count].find("'")
                            if first_sq_pos != -1:
                                second_sq_pos = hdr_list[count+tmp_count][first_sq_pos+1:].find("'")
                                if second_sq_pos != -1:
                                    new_entry = new_entry+hdr_list[count+tmp_count][first_sq_pos+1:second_sq_pos+1+first_sq_pos].rstrip()+"'"
                                else:
                                    # Raise exception here because there should be a second single quote after CONTINUE
                                    raise RuntimeError("There should be two single quotes after CONTINUE. Found the first, but not the second.")

                            else:
                                # Raise exception here because there should be a (first) single quote after CONTINUE
                                raise RuntimeError("There should be two single quotes after CONTINUE. Did not find any.")

                        tmp_count += 1

                    hdr_list_new.append(new_entry)

                else:
                    continue

        else:
            # Add END at the end of the header
            hdr_list_new.append(hdr_list[count])

    # Now we stitch together the CONTINUE information correctly,
    # with a "\n" at the end that we use as a separator later on
    # when we convert from a string to an astropy header.
    for count, item in enumerate(hdr_list_new):
        if len(item) > 80:
            this_entry = item[0:78]+"&'\n"
            rest       = "CONTINUE  '"+item[78:]
            while len(rest) > 80:
                this_entry = this_entry + rest[0:78] + "&'\n"
                rest       = "CONTINUE  '"+ rest[78:]
            this_entry = this_entry + rest
            hdr_list_new[count] = this_entry

    # Now we should have the correct list of strings. Since we can't convert a list to a
    # FITS header directly, we have to convert it to a string first, separated by "\n".
    hdr_str_new = '\n'.join([str(item) for item in hdr_list_new])

    # And finally we create the new corrected astropy FITS header from that string
    hdr_corr = fits.Header.fromstring(hdr_str_new, sep='\n')

    # Return the corrected header
    return hdr_corr


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
        hdr = sunkit_instruments.suvi.read_suvi(the_input, return_header_only=True)
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
