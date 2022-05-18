import astropy

from sunkit_instruments.suvi import io


def test_read_suvi_l1b_nc(L1B_NC):
    l1b_nc_header, l1b_nc_data, l1b_nc_dqf = io.read_suvi(L1B_NC)
    assert isinstance(l1b_nc_data, astropy.io.fits.header.header)
    assert l1b_nc_header.shape == (1280, 1280)
    assert l1b_nc_dqf.shape == (1280, 1280)


def test_read_suvi_l1b_fits(L1B_FITS):
    l1b_fits_header, l1b_fits_data, l1b_fits_dqf = io.read_suvi(L1B_FITS)
    assert isinstance(l1b_fits_header, astropy.io.fits.header.header)
    assert l1b_fits_data.shape == (1280, 1280)
    assert l1b_fits_dqf.shape == (1280, 1280)


def test_read_suvi_l2_composite(L2_COMPOSITE):
    l2_header, l2_data, _ = io.read_suvi(L2_COMPOSITE)
    assert isinstance(l2_header, astropy.io.fits.hdu.compressed.compimageheader)
    assert l2_data.shape == (1280, 1280)


def test_suvi_fix_l1b_header(L1B_FITS):
    header = io.fix_l1b_header(L1B_FITS)
    assert isinstance(header, astropy.io.fits.header.header)
