# Allow all possible file extensions for FITS and netCDF
FITS_FILE_EXTENSIONS = (
    ".fits",
    ".fts",
    ".fits.gz",
    "fts.gz",
)
NETCDF_FILE_EXTENSIONS = (
    ".nc",
    ".nc.gz",
    ".cdf",
    ".cdf.gz",
)
# Naming scheme for SUVI files
COMPOSITE_MATCHES = [
    "-l2-ci094",
    "-l2-ci131",
    "-l2-ci171",
    "-l2-ci195",
    "-l2-ci284",
    "-l2-ci304",
]
L1B_MATCHES = [
    "-L1b-Fe093",
    "-L1b-Fe131",
    "-L1b-Fe171",
    "-L1b-Fe195",
    "-L1b-Fe284",
    "-L1b-He303",
]
VALID_WAVELENGTH_CHANNELS = [94, 131, 171, 195, 284, 304]
VALID_SPACECRAFT = [16, 17]
FLIGHT_MODEL = {16: "FM1", 17: "FM2"}
# The setup for filter wheel 1 and 2 for the different exposure types
FILTER_SETUP = {
    94: {
        "long": {"FW1": "thin_zirconium", "FW2": "open"},
        "short": {"FW1": "thin_zirconium", "FW2": "open"},
        "short_flare": {"FW1": "thin_zirconium", "FW2": "thin_zirconium"},
    },
    131: {
        "long": {"FW1": "thin_zirconium", "FW2": "open"},
        "short": {"FW1": "thin_zirconium", "FW2": "open"},
        "short_flare": {"FW1": "thin_zirconium", "FW2": "thin_zirconium"},
    },
    171: {
        "long": {"FW1": "thin_aluminum", "FW2": "open"},
        "short_flare": {"FW1": "thin_aluminum", "FW2": "thin_aluminum"},
    },
    195: {
        "long": {"FW1": "thin_aluminum", "FW2": "open"},
        "short_flare": {"FW1": "thin_aluminum", "FW2": "thin_aluminum"},
    },
    284: {
        "long": {"FW1": "thin_aluminum", "FW2": "open"},
        "short_flare": {"FW1": "thin_aluminum", "FW2": "thin_aluminum"},
    },
    304: {
        "long": {"FW1": "thin_aluminum", "FW2": "open"},
        "short_flare": {"FW1": "thin_aluminum", "FW2": "thin_aluminum"},
    },
}
# This is how the global attributes of the netCDF file get
# mapped to the corresponding FITS header keywords.
TAG_MAPPING = {
    "instrument_id": "INST_ID",
    "platform_ID": "TELESCOP",
    "instrument_type": "INSTRUME",
    "project": "PROJECT",
    "institution": "ORIGIN",
    "production_site": "PRODSITE",
    "naming_authority": "NAMEAUTH",
    "production_environment": "PROD_ENV",
    "production_data_source": "DATA_SRC",
    "processing_level": "LEVEL",
    "algorithm_version": "CREATOR",
    "title": "TITLE",
    "keywords_vocabulary": "KEYVOCAB",
    "date_created": "DATE",
    "orbital_slot": "ORB_SLOT",
    "dataset_name": "FILENAME",
    "iso_series_metadata_id": "ISO_META",
    "id": "UUID",
    "LUT_Filenames": "LUT_NAME",
    "license": "LICENSE",
    "keywords": "KEYWORDS",
    "summary": "SUMMARY",
}
# Mapping for the FITS header keywords and their comments
TAG_COMMENT_MAPPING = {
    "SIMPLE": "file does conform to FITS standard",
    "BITPIX": "number of bits per data pixel",
    "NAXIS": "number of data axes",
    "NAXIS1": "length of data axis 1",
    "NAXIS2": "length of data axis 2",
    "EXTEND": "FITS dataset may contain extensions",
    "IMSENUMB": "[1] Image Serial Number",
    "CRPIX1": "[1] center of sun pixel in image along 1st axis",
    "CRPIX2": "[1] center of sun pixel in image along 2nd axis",
    "CDELT1": "[arcsec] 1st axis detector plate scale @ref pix",
    "CDELT2": "[arcsec] 2nd axis detector plate scale @ref pix",
    "DIAM_SUN": "[count] sun diameter in pixels",
    "CUNIT1": "1st axis detector plate scale units",
    "CUNIT2": "2nd axis detector plate scale units",
    "ORIENT": "orientation of image",
    "CROTA": "[degree] solar north pole angular offset",
    "SOLAR_B0": "[degree] solar equator angular offset",
    "PC1_1": "[1] 1st row, 1st col 2D transformation matrix",
    "PC1_2": "[1] 1st row, 2nd col 2D transformation matrix",
    "PC2_1": "[1] 2nd row, 1st col 2D transformation matrix",
    "PC2_2": "[1] 2nd row, 2nd col 2D transformation matrix",
    "CSYER1": "[arcsec] 1st axis systematic errors",
    "CSYER2": "[arcsec] 2nd axis systematic errors",
    "WCSNAME": "solar image coordinate system type",
    "CTYPE1": "1st axis coordinate system name",
    "CTYPE2": "2nd axis coordinate system name",
    "CRVAL1": "[degree] longitude of sun center for HPLN-TAN",
    "CRVAL2": "[degree] latitude of sun center for HPLT-TAN",
    "LONPOLE": "[degree] longitude of celestial north pole",
    "TIMESYS": "principal time system",
    "DATE-OBS": "sun observation start time on sat",
    "DATE-END": "sun observation end time on sat",
    "CMD_EXP": "[s] commanded imaging exposure time",
    "EXPTIME": "[s] actual imaging exposure time",
    "OBSGEO-X": "[m] observing platform ECEF X coordinate",
    "OBSGEO-Y": "[m] observing platform ECEF Y coordinate",
    "OBSGEO-Z": "[m] observing platform ECEF Z coordinate",
    "DSUN_OBS": "[m] distance to center of sun",
    "OBJECT": "object being viewed",
    "SCI_OBJ": "science objective of observation",
    "WAVEUNIT": "solar image wavelength units",
    "WAVELNTH": "[angstrom] solar image wavelength",
    "IMG_MIN": "[W m-2 sr-1] minimum radiance in image",
    "IMG_MAX": "[W m-2 sr-1] maximum radiance in image",
    "IMG_MEAN": "[W m-2 sr-1] mean radiance in image",
    "FILTER1": "forward filter setting mnemonic",
    "FILTER2": "aft filter setting mnemonic",
    "GOOD_PIX": "[count] number of good quality pixels in image",
    "FIX_PIX": "[count] number of corrected pixels in image",
    "SAT_PIX": "[count] number of saturated pixels in image",
    "MISS_PIX": "[count] number of missing pixels in image",
    "IMGTII": "[W m-2] total irradiance of image",
    "IMGTIR": "[W m-2 sr-1] total radiance of image",
    "IMG_SDEV": "[W m-2 sr-1] std dev of radiance in image",
    "EFF_AREA": "[m2] effective telescope area",
    "APSELPOS": "[1] aperture selector setting",
    "INSTRESP": "[count photon-1 cm-2] instrument response, used",
    "PHOT_ENG": "[J] photon energy, used in the calculation of r",
    "RSUN": "[count] solar angular radius in pixels",
    "HGLT_OBS": "[degree] Heliographic Stonyhurst Latitude of th",
    "HGLN_OBS": "[degree] Heliographic Stonyhurst Longitude of t",
    "HEEX_OBS": "[m] Heliocentric Earth Ecliptic X-axis coordina",
    "HEEY_OBS": "[m] Heliocentric Earth Ecliptic Y-axis coordina",
    "HEEZ_OBS": "[m] Heliocentric Earth Ecliptic Z-axis coordina",
    "FILTPOS1": "[1] forward filter wheel setting",
    "FILTPOS2": "[1] aft filter wheel setting",
    "YAW_FLIP": "[1] 0=upright 1=neither 2=inverted",
    "CCD_READ": "[1] CCD cnfg: 0=no cnfg 1=left amp 2=right amp",
    "ECLIPSE": "[1] sun obscured: 0=no eclipse 1=penumbra,prece",
    "CONTAMIN": "[angstrom] contamination thickness in angstroms",
    "CONT_FLG": "[1] contamination correction: 0=true 1=false",
    "DATE-BKE": "last contamination bake-out end time",
    "DER_SNR": "[W m-2 sr-1] CCD signal to noise ratio",
    "SAT_THR": "[W m-2 sr-1] CCD saturation point",
    "CCD_BIAS": "[count] CCD background electronic noise",
    "CCD_TMP1": "[degrees_C] sensor 1 camera temperature",
    "CCD_TMP2": "[degrees_C] sensor 2 camera temperature",
    "DATE-DFM": "median value dark frame time stamp",
    "NDFRAMES": "[count] number of source dark frames",
    "DATE-DF0": "1st observed dark frame time stamp",
    "DATE-DF1": "2nd observed dark frame time stamp",
    "DATE-DF2": "3rd observed dark frame time stamp",
    "DATE-DF3": "4th observed dark frame time stamp",
    "DATE-DF4": "5th observed dark frame time stamp",
    "DATE-DF5": "6th observed dark frame time stamp",
    "DATE-DF6": "7th observed dark frame time stamp",
    "DATE-DF7": "8th observed dark frame time stamp",
    "DATE-DF8": "9th observed dark frame time stamp",
    "DATE-DF9": "10th observed dark frame time stamp",
    "SOLCURR1": "[count] solar array current chan 1-4 in DN",
    "SOLCURR2": "[count] solar array current chan 5-8 in DN",
    "SOLCURR3": "[count] solar array current chan 9-12 in DN",
    "SOLCURR4": "[count] solar array current chan 13-16 in DN",
    "PCTL0ERR": "[percent] uncorrectable L0 error pct",
    "LONGSTRN": "The HEASARC Long String Convention may be used",
}
SOLAR_CLASSES = [
    ("unlabeled", 0),
    ("outer_space", 1),
    ("bright_region", 3),
    ("filament", 4),
    ("prominence", 5),
    ("coronal_hole", 6),
    ("quiet_sun", 7),
    ("limb", 8),
    ("flare", 9),
]
SOLAR_CLASS_NAME = {number: theme for theme, number in SOLAR_CLASSES}
SOLAR_COLORS = {
    "unlabeled": "white",
    "outer_space": "black",
    "bright_region": "#F0E442",
    "filament": "#D55E00",
    "prominence": "#E69F00",
    "coronal_hole": "#009E73",
    "quiet_sun": "#0072B2",
    "limb": "#56B4E9",
    "flare": "#CC79A7",
}
