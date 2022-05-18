import pytest

from sunkit_instruments.data.test import get_test_filepath


@pytest.fixture
def L1B_FITS():
    return get_test_filepath(
        "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz",
    )


@pytest.fixture
def L1B_NC():
    return get_test_filepath(
        "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc",
    )


@pytest.fixture
def L2_COMPOSITE():
    return get_test_filepath(
        "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits"
    )
