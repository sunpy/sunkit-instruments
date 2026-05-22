import numpy as np
import pytest

from astropy import units as u

from sunpy.map import Map

from sunkit_instruments import metis

_BASE_METIS_HEADER = {
    "INSTRUME": "METIS",
    "OBSRVTRY": "SOLAR ORBITER",
    "FILTER": "VL",
    "BTYPE": "VL total brightness",
    "LEVEL": "L2",
    "DATE-AVG": "2023-01-01T12:00:00",
    "CDELT1": 10.0,
    "CDELT2": 10.0,
    "CUNIT1": "arcsec",
    "CUNIT2": "arcsec",
    "CTYPE1": "HPLN-TAN",
    "CTYPE2": "HPLT-TAN",
    "NAXIS1": 1024,
    "NAXIS2": 1024,
    "CRPIX1": 512.0,
    "CRPIX2": 512.0,
    "CRVAL1": 0.0,
    "CRVAL2": 0.0,
    "RSUN_ARC": 960.0,
    "DSUN_OBS": 1.5e11,
    "HGLN_OBS": 0.0,
    "HGLT_OBS": 0.0,
    "INN_FOV": 1.6,
    "OUT_FOV": 2.9,
    "IO_XCEN": 512.0,
    "IO_YCEN": 512.0,
    "FS_XCEN": 512.0,
    "FS_YCEN": 512.0,
    "SUN_XCEN": 512.0,
    "SUN_YCEN": 512.0,
}

@pytest.fixture(scope="module")
def metis_map(metis_test_data):
    header = {**_BASE_METIS_HEADER}
    return Map(metis_test_data, header)

@pytest.fixture(scope="module")
def metis_test_data():
    """
    Generates synthetic 1024x1024 image data with a radial gradient.
    """
    size = 1024
    y, x = np.ogrid[-size / 2 : size / 2, -size / 2 : size / 2]
    r = np.sqrt(x**2 + y**2)
    # Exponential decay to simulate coronal brightness
    data = 1000 * np.exp(-r / 200) + np.random.normal(10, 2, (size, size))
    return data.astype(np.float32)

@pytest.fixture
def quality_matrix():
    """
    Generates a mock METIS Quality Matrix (QMatrix).
    Values: 1 = Good, 0 = Bad/Masked.
    """
    qmat = np.ones((1024, 1024), dtype=np.float32)
    # Simulate a block of bad pixels
    qmat[100:200, 100:200] = 0
    return qmat

"""Tests for METIS-specific data processing methods like qmatrix masking."""

def test_mask_bad_pix(metis_map, quality_matrix):
    """
    Test that mask_bad_pix correctly applies the Quality Matrix (QMatrix).
    """

    # Store the original value to check it remains unchanged later
    original_value = metis_map.data[512, 512]

    # Apply the quality matrix mask
    metis_map = metis.mask_bad_pix(metis_map = metis_map,qmat= quality_matrix)

    # 1. Pixels where quality_matrix was 0 should now be NaN
    assert np.isnan(metis_map.data[150, 150])

    # 2. Pixels where quality_matrix was 1 should still have their original value
    # and definitely should NOT be NaN
    assert not np.isnan(metis_map.data[512, 512])
    assert metis_map.data[512, 512] == original_value


def test_get_fov_rsun(metis_map):
    """
    Test the get_fov_rsun method returns correct solar radii values.
    """

    # Metis specific method to get FOV in R_sun
    rmin, rmax, board = metis_map._get_fov_rsun()

    assert isinstance(rmin, u.Quantity)
    assert rmin > 0
    assert rmax > rmin
    # With INN_FOV=1.6 deg and RSUN_ARC=960", rmin should be ~6.0 R_sun
    # (1.6 * 3600 / 960 = 6.0)
    assert pytest.approx(rmin, rel=1e-2) == 6.0

def test_mask_bad_pix_invalid_input(metis_map):
    """Verify that mask_bad_pix raises errors on wrong dimensions."""
    wrong_shape_qmat = np.ones((512, 512))

    with pytest.raises(ValueError, match="shape"):
        metis.mask_bad_pix(metis_map,wrong_shape_qmat)


def test_mask_occs_logic(metis_map):
    """
    Ensure the mask_occs method successfully sets internal occultor
    pixels to NaN.
    """
    metis_map = metis.mask_occs(metis_map)
    assert np.isnan(metis_map.data[512, 512])

@pytest.mark.parametrize("mask_val", [np.nan, -999])
def test_mask_occs_various_values(metis_map, mask_val):
    metis_map = metis.mask_occs(metis_map,mask_val)
    if np.isnan(mask_val):
        assert np.sum(np.isnan(metis_map.data)) > 0
    else:
        assert np.sum(metis_map.data == mask_val) > 0
