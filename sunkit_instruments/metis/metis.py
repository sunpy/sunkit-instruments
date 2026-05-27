
import copy

import numpy as np

from astropy import units as u
from astropy.visualization import ImageNormalize

from sunpy.map.sources.solo import METISMap
from sunpy.util.exceptions import warn_user


def mask_occs(metis_map: METISMap, mask_val : np.ndarray = np.nan) -> METISMap:
    """
    Mask the data in regions obscured by internal and external occulters.

    This method modifies the data array in-place, setting pixels outside
    the field of view to the specified mask value.

    Parameters
    ----------
    mask_val : `float`, optional
        The value to assign to masked pixels (outside the field of view).
        Default is ``np.nan``.

    Notes
    -----
    This method directly modifies ``metis_map.data`` in-place.

    Warnings
    --------
    If CDELT1 and CDELT2 are not equal (non-square pixels), this method
    will warn and exit without performing masking, as the circular mask
    calculation assumes square pixels.

    """

    metis_map = copy.deepcopy(metis_map)

    if metis_map.scale[0] != metis_map.scale[1]:
        warn_user(
            "Warning: CDELT1 != CDELT2 for {fname}. Exiting mask_occs method...".format(fname=metis_map.meta["filename"])
        )
        return


    # Calculate occulter radii in pixels
    inn_fov = (metis_map.meta["inn_fov"] * u.deg / metis_map.scale[0]).decompose()  # in pix
    out_fov = (metis_map.meta["out_fov"] * u.deg / metis_map.scale[1]).decompose()  # in pix

    # Create coordinate grids
    x = np.arange(0, metis_map.meta["naxis1"], 1)
    y = np.arange(0, metis_map.meta["naxis2"], 1)
    xx, yy = np.meshgrid(x, y, sparse=True)

    # Calculate distance from internal occulter center
    in_xcen = (metis_map.meta["io_xcen"] - 1) * u.pix
    in_ycen = (metis_map.meta["io_ycen"] - 1) * u.pix
    dist_inncen = np.sqrt((xx * u.pix - in_xcen) ** 2 + (yy * u.pix - in_ycen) ** 2)

    # Calculate distance from external occulter/field stop center
    # NOTE: Workaround for DR1 data where fs_*cen keywords are not correctly defined.
    # In DR1, fs_xcen and fs_ycen are not available and as consequence the corresponding keywords in a FITS file are assigned to crpix1 and crpix2, respectively.
    # When this occurs, sun_xcen and sun_ycen should be used as approximate coordinates of the field stop center.
    # This problem is solved in DR2.
    if metis_map.meta["fs_xcen"] == metis_map.meta["crpix1"] and metis_map.meta["fs_ycen"] == metis_map.meta["crpix2"]:
        # DR1 workaround: use sun center instead; For the DR1 data fs_*cen keywords are not defined correctly in a FITS file header
        out_xcen = (metis_map.meta["sun_xcen"] - 1) * u.pix
        out_ycen = (metis_map.meta["sun_ycen"] - 1) * u.pix
    else:
        # Normal case: use field stop center
        out_xcen = (metis_map.meta["fs_xcen"] - 1) * u.pix
        out_ycen = (metis_map.meta["fs_ycen"] - 1) * u.pix

    dist_outcen = np.sqrt((xx * u.pix - out_xcen) ** 2 + (yy * u.pix - out_ycen) ** 2)
    # Apply masks
    metis_map.data[dist_inncen <= inn_fov] = mask_val
    metis_map.data[dist_outcen >= out_fov] = mask_val
    
    return metis_map

def mask_bad_pix(metis_map:METISMap, qmat:np.ndarray , mask_val=np.nan ) -> METISMap:
    """
    Mask bad-quality pixels in the Metis image.

    This method modifies the data array, setting bad-quality
    pixels to the specified mask value based on the provided quality matrix.

    Parameters
    ----------
    qmat : `numpy.ndarray`
        Pixel quality matrix with the same shape as the image data.
        Expected values:
        - ``1`` : linear range (good-quality pixels)
        - ``0`` : close to 0 counts or close to saturation (bad-quality)
        - ``np.nan`` : exactly 0 count or saturated pixels (bad-quality)
    mask_val : `float`, optional
        The value to assign to masked bad pixels. Default is ``np.nan``.

    Raises
    ------
    ValueError
        If the quality matrix shape does not match the data shape.
    TypeError
        If qmat is not a numpy array.

    """
    metis_map = copy.deepcopy(metis_map)
    # Validate input type
    if not isinstance(qmat, np.ndarray):
        raise TypeError(f"qmat must be a numpy.ndarray, got {type(qmat).__name__}")
    # Validate shape compatibility
    if qmat.shape != metis_map.data.shape:
        raise ValueError(
            f"Pixel quality matrix shape {qmat.shape} does not match "
            f"METISMap data shape {metis_map.data.shape}. Cannot apply mask."
        )
    # Create mask: keep only pixels with value 1 (good quality)
    qmat_mask = qmat == 1
    metis_map.data[~qmat_mask] = mask_val

    return metis_map
