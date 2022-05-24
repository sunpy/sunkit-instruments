import os
import gzip
import logging
import tempfile
from pathlib import Path

import h5py
import numpy

import sunpy.map
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from sunpy.util.exceptions import warn_user

from sunkit_instruments.suvi._variables import (
    COMPOSITE_MATCHES,
    FITS_FILE_EXTENSIONS,
    L1B_MATCHES,
    NETCDF_FILE_EXTENSIONS,
    TAG_COMMENT_MAPPING,
    TAG_MAPPING,
)

__all__ = ["read_suvi", "files_to_map"]


def _fix_l1b_header(filename):
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
    filename: `str`
        Filename of the L1b file with the corrupt FITS header.

    Returns
    -------
    `astropy.io.fits.header.Header`
        Corrected FITS header.
    """
    try:
        # First try it with the astropy .to_string() method, as this is the easiest.
        hdr = fits.getheader(filename)
        hdr_str = hdr.tostring()
    except Exception:
        # Read the file manually as bytes until we hit a UnicodeDecodeError, i.e.
        # until we reach the data part. Since astropy version 4.2.1, we can't use
        # the .to_string() method anymore because of FITS header consistency checks
        # that cannot be overridden, and they won't fix it unfortunately. If the
        # input file is a .gz file, we need to unpack it first to the tmp directory.
        temp_dir = tempfile.gettempdir()
        name = Path(filename).name
        is_gz_file = False
        if name.endswith(".gz"):
            is_gz_file = True
            with gzip.open(filename, "r") as gfile:
                filename = str(Path(temp_dir) / name[:-3])
                with open(filename, "wb") as file_out:
                    file_out.write(gfile.read())
        hdr_str = ""
        with open(filename, "rb") as file:
            counter = 1
            while True:
                try:
                    this_line = file.read(counter)
                    this_str = this_line.decode("utf-8")
                    hdr_str += this_str
                    counter += 1
                except UnicodeDecodeError:
                    break
        if is_gz_file:
            os.remove(filename)
    # Make a list of strings with a length of 80
    hdr_list = [hdr_str[i : i + 80] for i in range(0, len(hdr_str), 80)]
    # Remove all the empty entries
    while " " * 80 in hdr_list:
        hdr_list.remove(" " * 80)
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
                                    raise RuntimeError(
                                        "There should be two single quotes after CONTINUE. Found the first, but not the second."
                                    )
                            else:
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
    hdr_corr = fits.Header.fromstring(hdr_str_new, sep="\n")
    return hdr_corr


def _read_fits(filename):
    """
    Read a FITS file and return the header, data and dqf.
    """
    if any(fn in os.path.basename(filename) for fn in COMPOSITE_MATCHES):
        with fits.open(filename) as hdu:
            data, header = hdu[1].data, hdu[1].header
            dqf = None
    elif any(fn in os.path.basename(filename) for fn in L1B_MATCHES):
        with fits.open(filename) as hdu:
            data, header, dqf = hdu[0].data, _fix_l1b_header(filename), hdu[1].data
    else:
        raise ValueError(
            f"File {filename} does not look like a SUVI L1b FITS file or L2 HDR composite."
        )
    return header, data, dqf


def _make_cdf_header(header_info):
    header_info_copy = header_info.copy()
    # Discard everything where the key name is longer than 8 characters,
    # plus specific entries we have to deal with manually.
    for key, value in header_info.items():
        if len(key) > 8:
            del header_info_copy[key]
        elif key in ["RAD", "DQF", "NAXIS1", "NAXIS2"]:
            del header_info_copy[key]
    for key, value in header_info_copy.items():
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
                    header_info_copy[key] = numpy.longlong(value)
                elif value.dtype in [numpy.float16, numpy.float32, numpy.float64]:
                    header_info_copy[key] = numpy.float64(value)
            else:
                if value.dtype == "|S1":
                    # Byte string to actual string, and removing weird characters
                    header_info_copy[key] = (
                        value.tobytes().decode("utf-8").rstrip("\x00")
                    )
    # Now deal with the dates (float in the netCDF). Transform to readable string,
    # ignore bakeout date because it is always -999.
    for key, value in header_info_copy.items():
        if key.startswith("DATE") and key != "DATE-BKE":
            # Explanation for the odd time creation: the SUVI files say they use the
            # the J2000 epoch, but they do not: the reference time is 2000-01-01 at
            # 12:00:00 *UTC*, whereas the reference time for J2000 is in *TT*. So in
            # order to get the time right, we need to define it in TT, but add the
            # offset of 69.184 seconds between UTC and TT.
            the_readable_date = (
                Time("2000-01-01T12:01:09.184", scale="tt") + value * u.s
            )
            header_info_copy[key] = the_readable_date.utc.value
    # Add NAXIS1 and NAXIS2 manually, because they are odd coming from the netCDF
    header_info_copy["NAXIS1"] = None
    header_info_copy["NAXIS2"] = None
    # Same for BLANK, BSCALE, and BZERO
    header_info_copy["BLANK"] = None
    header_info_copy["BSCALE"] = None
    header_info_copy["BZERO"] = None
    header_info_copy["BUNIT"] = None
    header = fits.Header.fromkeys(header_info_copy.keys())
    for keyword in header:
        header[keyword] = header_info_copy[keyword]
    # Add fits header comments for known keywords as defined above
    for keyword in header:
        if keyword in TAG_COMMENT_MAPPING:
            header.set(keyword, header[keyword], TAG_COMMENT_MAPPING[keyword])
    # Add EXTEND, EXTVER, EXTNAME, and LONGSTR
    header.append(("EXTEND", True, "FITS dataet may contain extensions"))
    header.append(("EXTVER", 1, ""))
    header.append(("EXTNAME", "DATA", ""))
    header.append(
        ("LONGSTRN", "OGIP 1.0", "The HEASARC Long String Convention may be used")
    )
    return header


def _read_netCDF(filename):
    """
    Read a CDF file and return the header, data and dqf.
    """
    if any(fn in os.path.basename(filename) for fn in L1B_MATCHES):
        with h5py.File(filename, "r") as afile:
            data = afile["RAD"][:]

            blank = afile["RAD"].attrs["_FillValue"][0]
            bzero = afile["RAD"].attrs["add_offset"][0]
            bscale = afile["RAD"].attrs["scale_factor"][0]
            bunit = afile["RAD"].attrs["units"].tobytes().decode("utf-8").rstrip("\x00")

            data = data * bscale + bzero
            dqf = afile["DQF"][:]

            header_info = dict((key, afile[key][...]) for key in afile.keys())
            header = _make_cdf_header(header_info)
            # Deal with this here as we require the file.
            for att, val in afile.attrs.items():
                if att in TAG_MAPPING:
                    header[TAG_MAPPING[att]] = (
                        val.tobytes().decode("utf-8").rstrip("\x00")
                    )
            header["NAXIS1"] = data.shape[0]
            header["NAXIS2"] = data.shape[1]
            header["BLANK"] = blank
            header["BSCALE"] = bscale
            header["BZERO"] = bzero
            header["BUNIT"] = bunit
    else:
        raise ValueError(f"File {filename} does not look like a SUVI L1b netCDF file.")
    return header, data, dqf


def read_suvi(filename):
    """
    Read a SUVI L1b FITS or netCDF file or a L2 HDR composite FITS file.

    Returns header, data and the data quality flag array (DQF) for L1b files.

    For SUVI L1b FITS files, the broken FITS header is fixed automatically
    (broken because of the wrong implementation of the CONTINUE convention).

    This read function is intended to provide a consistent file interface
    for FITS and netCDF, L1b and L2.

    .. note::
        The type of file is determined by pattern matching in the
        filenames, e.g. "-L1b-Fe171" for a 171 L1b file and "-l2-ci171"
        for a 171 L2 HDR composite. If those patterns are not found
        in the filename, the files will not be recognized.

    .. note::
        If ``filename`` is an L1b netCDF file, the information from
        the netCDF file is transformed into a FITS header.

    Parameters
    ----------
    filename : `str`
        File to read.

    Returns
    -------
    `astropy.io.fits.header.Header`, `~numpy.ndarray`, `~numpy.ndarray`
        Header, data, and data quality flags.
    """
    if filename.lower().endswith(FITS_FILE_EXTENSIONS):
        header, data, dqf = _read_fits(filename)
    elif filename.lower().endswith(NETCDF_FILE_EXTENSIONS):
        header, data, dqf = _read_netCDF(filename)
    else:
        raise ValueError(
            f"File {filename} does not look like a valid FITS or netCDF file."
        )
    return header, data, dqf


def files_to_map(
    files,
    despike_l1b=False,
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
    despike_l1b: `bool`, optional. Default: False.
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
        A map (sequence) of the SUVI data, or `None` if no
        data was found matching the given criteria.
    """
    # Avoid circular imports
    from sunkit_instruments.suvi.suvi import despike_l1b_array

    if isinstance(files, str):
        files = [files]
    files = sorted(files)
    if any(fn in os.path.basename(files[0]) for fn in COMPOSITE_MATCHES):
        composites = True
    elif any(fn in os.path.basename(files[0]) for fn in L1B_MATCHES):
        composites = False
    else:
        raise ValueError(
            f"First file {files[0]} does not look like a SUVI L1b file or L2 HDR composite."
        )

    datas = []
    headers = []
    for afile in files:
        logging.debug(f"Reading {afile}")
        if composites:
            if any(fn in os.path.basename(afile) for fn in COMPOSITE_MATCHES):
                header, data, _ = read_suvi(afile)
                datas.append(data)
                headers.append(header)
            else:
                warn_user(
                    f"File {afile} does not look like a SUVI L2 HDR composite. Skipping."
                )
        else:
            if any(fn in os.path.basename(afile) for fn in L1B_MATCHES):
                header, data, dqf_mask = read_suvi(afile)
                if despike_l1b:
                    data = despike_l1b_array(data, dqf_mask)
                if only_long_exposures:
                    if "long_exposure" in header["SCI_OBJ"]:
                        datas.append(data)
                        headers.append(header)
                elif only_short_exposures:
                    if "short_exposure" in header["SCI_OBJ"]:
                        datas.append(data)
                        headers.append(header)
                elif only_short_flare_exposures:
                    if "short_flare_exposure" in header["SCI_OBJ"]:
                        datas.append(data)
                        headers.append(header)
                else:
                    datas.append(data)
                    headers.append(header)
            else:
                warn_user(f"File {afile} does not look like a SUVI L1b file. Skipping.")
    if len(datas) == 1:
        return sunpy.map.Map(datas[0], headers[0])
    elif len(datas) > 1:
        return sunpy.map.Map(list(zip(datas, headers)), sequence=True)
    else:
        warn_user("List of data/headers is empty.")
