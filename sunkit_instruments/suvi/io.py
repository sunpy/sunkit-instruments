import os
import gzip
import tempfile

import h5py
import matplotlib.pyplot as plt
import numpy
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

import sunpy.map
from astropy import units as u
from astropy.io import fits
from astropy.io.fits.verify import VerifyError
from astropy.time import Time
from sunpy.util.exceptions import warn_user

from sunkit_instruments.suvi.io import read_suvi

__all__ = ["fix_L1b_header", "read_suvi"]


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
        if split_filename[1] == ".gz":
            is_gz_file = True
            with gzip.open(input_filename, "r") as f_in, open(
                temp_dir + split_filename[0], "wb"
            ) as f_out:
                f_out.write(f_in.read())
            file_to_open = temp_dir + split_filename[0]
        else:
            is_gz_file = False
            file_to_open = input_filename

        hdr_str = ""
        with open(file_to_open, "rb") as in_file:
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
    hdr_list = [hdr_str[i : i + 80] for i in range(0, len(hdr_str), 80)]

    # Remove all the empty entries
    while " " * 80 in hdr_list:
        hdr_list.remove(" " * 80)

    # Make a new string list where we put all the information together correctly
    hdr_list_new = []
    for count, item in enumerate(hdr_list):
        if count <= len(hdr_list) - 2:
            if (
                hdr_list[count][0:8] != "CONTINUE"
                and hdr_list[count + 1][0:8] != "CONTINUE"
            ):
                hdr_list_new.append(hdr_list[count])
            else:
                if (
                    hdr_list[count][0:8] != "CONTINUE"
                    and hdr_list[count + 1][0:8] == "CONTINUE"
                ):
                    ampersand_pos = hdr_list[count].find("&")
                    if ampersand_pos != -1:
                        new_entry = hdr_list[count][0:ampersand_pos]
                    else:
                        # Raise exception here because there should be an ampersand at the end of a CONTINUE'd keyword
                        raise RuntimeError(
                            "There should be an ampersand at the end of a CONTINUE'd keyword."
                        )

                    tmp_count = 1
                    while hdr_list[count + tmp_count][0:8] == "CONTINUE":
                        ampersand_pos = hdr_list[count + tmp_count].find("&")
                        if ampersand_pos != -1:
                            first_sq_pos = hdr_list[count + tmp_count].find("'")
                            if first_sq_pos != -1:
                                new_entry = (
                                    new_entry
                                    + hdr_list[count + tmp_count][
                                        first_sq_pos + 1 : ampersand_pos
                                    ]
                                )
                            else:
                                # Raise exception here because there should be a single quote after CONTINUE
                                raise RuntimeError(
                                    "There should be two single quotes after CONTINUE. Did not find any."
                                )

                        else:
                            # If there is no ampersand at the end anymore, it means the entry ends here.
                            # Read from the first to the second single quote in this case.
                            first_sq_pos = hdr_list[count + tmp_count].find("'")
                            if first_sq_pos != -1:
                                second_sq_pos = hdr_list[count + tmp_count][
                                    first_sq_pos + 1 :
                                ].find("'")
                                if second_sq_pos != -1:
                                    new_entry = (
                                        new_entry
                                        + hdr_list[count + tmp_count][
                                            first_sq_pos
                                            + 1 : second_sq_pos
                                            + 1
                                            + first_sq_pos
                                        ].rstrip()
                                        + "'"
                                    )
                                else:
                                    # Raise exception here because there should be a second single quote after CONTINUE
                                    raise RuntimeError(
                                        "There should be two single quotes after CONTINUE. Found the first, but not the second."
                                    )

                            else:
                                # Raise exception here because there should be a (first) single quote after CONTINUE
                                raise RuntimeError(
                                    "There should be two single quotes after CONTINUE. Did not find any."
                                )

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
            this_entry = item[0:78] + "&'\n"
            rest = "CONTINUE  '" + item[78:]
            while len(rest) > 80:
                this_entry = this_entry + rest[0:78] + "&'\n"
                rest = "CONTINUE  '" + rest[78:]
            this_entry = this_entry + rest
            hdr_list_new[count] = this_entry

    # Now we should have the correct list of strings. Since we can't convert a list to a
    # FITS header directly, we have to convert it to a string first, separated by "\n".
    hdr_str_new = "\n".join([str(item) for item in hdr_list_new])

    # And finally we create the new corrected astropy FITS header from that string
    hdr_corr = fits.Header.fromstring(hdr_str_new, sep="\n")

    # Return the corrected header
    return hdr_corr


# read_fits helper function
def _read_fits(input_filename, return_DQF=False, return_header_only=False):
    data = None
    dqf = None

    # Based on the filename, determine if we are dealing
    # with L1b files or HDR composites (or neither).
    if any(fn in os.path.basename(input_filename) for fn in COMPOSITE_MATCHES):
        if return_header_only:
            header = fits.getheader(input_filename, 1)
        else:
            hdu = fits.open(input_filename)
            data, header = hdu[1].data, hdu[1].header
            hdu.close()
    elif any(fn in os.path.basename(input_filename) for fn in L1B_MATCHES):
        hdu = fits.open(input_filename)
        if return_header_only:
            header = fix_L1b_header(input_filename)
        else:
            hdu = fits.open(input_filename)
            data, header = hdu[0].data, fix_L1b_header(input_filename)
        if return_DQF:
            dqf = hdu[1].data
        hdu.close()
    else:
        raise ValueError(
            "File "
            + input_filename
            + " does not look like a SUVI L1b FITS file or L2 HDR composite."
        )

    return {"header": header, "data": data, "dqf": dqf}


# read_netCDF helper function
def _read_netCDF(input_filename, return_DQF=False, return_header_only=False):
    data = None
    dqf = None
    if any(fn in os.path.basename(input_filename) for fn in L1B_MATCHES):
        tmp_file = h5py.File(input_filename, "r")
        # Get the data first
        data = tmp_file["RAD"][:]
        # Get BLANK, BZERO, BSCALE, and BUNIT from the RAD attributes
        blank = tmp_file["RAD"].attrs["_FillValue"][0]
        bzero = tmp_file["RAD"].attrs["add_offset"][0]
        bscale = tmp_file["RAD"].attrs["scale_factor"][0]
        bunit = tmp_file["RAD"].attrs["units"].tobytes().decode("utf-8").rstrip("\x00")
        # Multiply the data accordingly
        data = data * bscale + bzero
        # Get the DQF if requested
        if return_DQF:
            dqf = tmp_file["DQF"][:]
        # Now deal with the header. Create a dictionary from the
        # the netCDF keys, and a copy of it.
        d = dict((key, tmp_file[key][...]) for key in tmp_file.keys())
        tmp_d = d.copy()
        # Discard everything where the key name is longer than 8 characters,
        # plus specific entries we have to deal with manually.
        for key, value in d.items():
            if len(key) > 8:
                del tmp_d[key]
            elif key in ["RAD", "DQF", "NAXIS1", "NAXIS2"]:
                del tmp_d[key]
        for key, value in tmp_d.items():
            if isinstance(value, numpy.ndarray):
                # We only want single values for the header, no arrays of length 1.
                # We convert everything that looks like an integer to a long,
                # everything that looks like a float to float64, and byte strings
                # to actual strings.
                if value.ndim == 0:
                    if value.dtype in [
                        numpy.int8,
                        numpy.int16,
                        numpy.int32,
                        numpy.int64,
                        numpy.uint8,
                        numpy.uint16,
                        numpy.uint32,
                        numpy.uint64,
                    ]:
                        tmp_d[key] = numpy.longlong(value)
                    elif value.dtype in [numpy.float16, numpy.float32, numpy.float64]:
                        tmp_d[key] = numpy.float64(value)
                else:
                    if value.dtype == "|S1":
                        # Byte string to actual string, and removing weird characters
                        tmp_d[key] = value.tobytes().decode("utf-8").rstrip("\x00")
        # Now deal with the dates (float in the netCDF). Transform to readable string,
        # ignore bakeout date because it is always -999.
        for key, value in tmp_d.items():
            if key.startswith("DATE") and key != "DATE-BKE":
                # Explanation for the odd time creation: the SUVI files say they use the
                # the J2000 epoch, but they do not: the reference time is 2000-01-01 at
                # 12:00:00 *UTC*, whereas the reference time for J2000 is in *TT*. So in
                # order to get the time right, we need to define it in TT, but add the
                # offset of 69.184 seconds between UTC and TT.
                the_readable_date = (
                    Time("2000-01-01T12:01:09.184", scale="tt") + value * u.s
                )
                tmp_d[key] = the_readable_date.utc.value
        # Add NAXIS1 and NAXIS2 manually, because they are odd coming from the netCDF
        tmp_d["NAXIS1"] = None
        tmp_d["NAXIS2"] = None
        # Same for BLANK, BSCALE, and BZERO
        tmp_d["BLANK"] = None
        tmp_d["BSCALE"] = None
        tmp_d["BZERO"] = None
        tmp_d["BUNIT"] = None
        header = fits.Header.fromkeys(tmp_d.keys())
        for keyword in header:
            header[keyword] = tmp_d[keyword]
        # Add all the info from the global attributes
        for att, val in tmp_file.attrs.items():
            if att in TAG_MAPPING:
                header[TAG_MAPPING[att]] = val.tobytes().decode("utf-8").rstrip("\x00")
        tmp_file.close()
        # And some manual info
        header["NAXIS1"] = data.shape[0]
        header["NAXIS2"] = data.shape[1]
        header["BLANK"] = blank
        header["BSCALE"] = bscale
        header["BZERO"] = bzero
        header["BUNIT"] = bunit
        # Add fits header comments for known keywords as defined above
        for keyword in header:
            if keyword in TAG_COMMENT_MAPPING:
                header.set(keyword, header[keyword], TAG_COMMENT_MAPPING[keyword])
        # Add EXTEND, EXTVER, EXTNAME, and LONGSTR
        header.append(("EXTEND", True, "FITS dataset may contain extensions"))
        header.append(("EXTVER", 1, ""))
        header.append(("EXTNAME", "DATA", ""))
        header.append(
            ("LONGSTRN", "OGIP 1.0", "The HEASARC Long String Convention may be used")
        )
    else:
        raise ValueError(
            "File " + input_filename + " does not look like a SUVI L1b netCDF file."
        )

    return {"header": header, "data": data, "dqf": dqf}


def read_suvi(input_filename, return_DQF=False, return_header_only=False):
    """
    Read a SUVI L1b FITS or netCDF file or a L2 HDR composite FITS file.
    Return data and header, optionally the data quality flag array (DQF)
    for L1b files. For SUVI L1b FITS files, the broken FITS header is
    fixed automatically (broken because of the wrong implementation of
    the CONTINUE convention). This read function is intented to provide
    a consistent file interface for FITS and netCDF, L1b and L2.

    .. note::
        The type of file is determined by pattern matching in the
        filenames, e.g. "-L1b-Fe171" for a 171 L1b file and "-l2-ci171"
        for a 171 L2 HDR composite. If those patterns are not found
        in the filename, the files will not be recognized.

    .. note::
        If input_filename is an L1b netCDF file, the information from
        the netCDF file is transformed into a FITS header.

    Parameters
    ----------
    input_filename: `str`
        File to read.

    return_DQF: `bool`, optional. Default: False.
        If True, returns the data quality flag array (DQF, only for L1b files), otherwise `None` for dqf.

    return_header_only: `bool`, optional. Default: False.
        If True, does not read the data array and returns `None` for data (the DQF can still be requested).

    Returns
    -------
    header, data, dqf: `astropy.io.fits.header.Header`, `~numpy.ndarray` (or `None` if return_header_only was set), `~numpy.ndarray` (or `None` if return_DQF was not set)
        header, data, and data quality flags.
    """

    # If it is a fits file...
    if input_filename.lower().endswith(FITS_FILE_EXTENSIONS):
        file_info = _read_fits(
            input_filename, return_DQF=return_DQF, return_header_only=return_header_only
        )

    # If it is a netCDF file...
    elif input_filename.lower().endswith(NETCDF_FILE_EXTENSIONS):
        file_info = _read_netCDF(
            input_filename, return_DQF=return_DQF, return_header_only=return_header_only
        )

    else:
        raise ValueError(
            "File "
            + input_filename
            + " does not look like a valid FITS or netCDF file."
        )

    return file_info["header"], file_info["data"], file_info["dqf"]


def files_to_map(
    files,
    sort_files=True,
    verbose=False,
    despike_L1b=False,
    only_long_exposures=False,
    only_short_exposures=False,
    only_short_flare_exposures=False,
):
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
    composite_matches = [
        "-l2-ci094",
        "-l2-ci131",
        "-l2-ci171",
        "-l2-ci195",
        "-l2-ci284",
        "-l2-ci304",
    ]
    L1b_matches = [
        "-L1b-Fe093",
        "-L1b-Fe131",
        "-L1b-Fe171",
        "-L1b-Fe195",
        "-L1b-Fe284",
        "-L1b-He303",
    ]
    if any(fn in os.path.basename(files[0]) for fn in composite_matches):
        composites = True
    elif any(fn in os.path.basename(files[0]) for fn in L1b_matches):
        composites = False
    else:
        raise ValueError(
            "First file "
            + files[0]
            + " does not look like a SUVI L1b file or L2 HDR composite."
        )

    datas = []
    headers = []
    for tmp_file in files:
        if verbose:
            print("Reading", tmp_file)
        # Test for L1b or composite based on the filename
        if composites:
            if any(fn in os.path.basename(tmp_file) for fn in composite_matches):
                tmp_header, tmp_data, _ = sunkit_instruments.suvi.read_suvi(tmp_file)
                datas.append(tmp_data)
                headers.append(tmp_header)
            else:
                warn_user(
                    "File "
                    + tmp_file
                    + " does not look like a SUVI L2 HDR composite. Skipping."
                )
        else:
            if any(fn in os.path.basename(tmp_file) for fn in L1b_matches):
                if despike_L1b:
                    tmp_header, tmp_data, dqf_mask = sunkit_instruments.suvi.read_suvi(
                        tmp_file, return_DQF=True
                    )
                    tmp_data = despike_L1b_image((tmp_data, dqf_mask))
                else:
                    tmp_header, tmp_data, _ = sunkit_instruments.suvi.read_suvi(
                        tmp_file
                    )
                if only_long_exposures:
                    if "long_exposure" in tmp_header["SCI_OBJ"]:
                        datas.append(tmp_data)
                        headers.append(tmp_header)
                elif only_short_exposures:
                    if "short_exposure" in tmp_header["SCI_OBJ"]:
                        datas.append(tmp_data)
                        headers.append(tmp_header)
                elif only_short_flare_exposures:
                    if "short_flare_exposure" in tmp_header["SCI_OBJ"]:
                        datas.append(tmp_data)
                        headers.append(tmp_header)
                else:
                    datas.append(tmp_data)
                    headers.append(tmp_header)
            else:
                warn_user(
                    "File "
                    + tmp_file
                    + " does not look like a SUVI L1b file. Skipping."
                )

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


def plot_thematic_map(input_filename, timestamp=True, legend=True, figsize=(10, 10)):
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

    if "-l2-thmap" in os.path.basename(input_filename):
        hdu = fits.open(input_filename)
        thmap_data = hdu[0].data
        if timestamp:
            time_st = hdu[0].header["DATE-OBS"][0:19]
        hdu.close()

        colortable = [
            SOLAR_COLORS[SOLAR_CLASS_NAME[i]] if i in SOLAR_CLASS_NAME else "black"
            for i in range(max(list(SOLAR_CLASS_NAME.keys())) + 1)
        ]
        cmap = ListedColormap(colortable)

        fig, ax = plt.subplots(constrained_layout=True, figsize=figsize)
        ax.imshow(
            thmap_data,
            origin="lower",
            cmap=cmap,
            vmin=-1,
            vmax=len(colortable),
            interpolation="none",
        )
        ax.set_axis_off()
        if timestamp:
            ax.text(10, 1245, time_st, fontsize=14, color="white")
        if legend:
            legend_elements = [
                Patch(facecolor=color, edgecolor="black", label=label.replace("_", " "))
                for label, color in SOLAR_COLORS.items()
            ]
            ax.legend(
                handles=legend_elements,
                loc="upper right",
                ncol=3,
                fancybox=True,
                shadow=False,
            )
        fig.show()
    else:
        raise ValueError(
            "File " + input_filename + " does not look like a SUVI L2 Thematic Map."
        )
