import shutil
import os.path

import pytest
import requests

import astropy

from sunkit_instruments import suvi
from sunkit_instruments.data.test import rootdir


@pytest.mark.remote_data
def test_io_download_suvi_testdata(tmp_path):
    try:
        f1 = requests.get("https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/"+
                          "suvi-l1b-fe171/2021/12/31/OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")
        open(os.path.join(str(tmp_path), "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz"), 'wb').write(f1.content)
    finally:
        f1.close()
    outfile1 = tmp_path / "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz"

    try:
        f2 = requests.get("https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/"+
                          "suvi-l2-ci171/2021/12/31/dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits")
        open(os.path.join(str(tmp_path), "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits"), 'wb').write(f2.content)
    finally:
        f2.close()
    outfile2 = tmp_path / "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits"

    try:
        f3 = requests.get("https://noaa-goes16.s3.amazonaws.com/SUVI-L1b-Fe171/2021/365/00/"+
                          "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc")
        open(os.path.join(str(tmp_path), "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc"), 'wb').write(f3.content)
    finally:
        f3.close()
    outfile3 = tmp_path / "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc"

    dst1 = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")
    dst2 = os.path.join(rootdir, "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits")
    dst3 = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc")

    shutil.copyfile(str(outfile1), dst1)
    shutil.copyfile(str(outfile2), dst2)
    shutil.copyfile(str(outfile3), dst3)

    assert os.path.exists(dst1)
    assert os.path.exists(dst2)
    assert os.path.exists(dst3)


def test_io_read_suvi_files():
    L1b_nc = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc")
    L1b_fits = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")
    L2_composite = os.path.join(rootdir, "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits")

    L1b_nc_header, L1b_nc_data, L1b_nc_dqf = suvi.read_suvi(L1b_nc, return_DQF=True)
    L1b_fits_header, L1b_fits_data, L1b_fits_dqf = suvi.read_suvi(L1b_fits, return_DQF=True)
    L2_header, L2_data = suvi.read_suvi(L2_composite)

    assert isinstance(L1b_nc_header, astropy.io.fits.header.Header)
    assert isinstance(L1b_fits_header, astropy.io.fits.header.Header)
    assert isinstance(L2_header, astropy.io.fits.hdu.compressed.CompImageHeader)

    assert L1b_nc_data.shape == (1280, 1280)
    assert L1b_fits_data.shape == (1280, 1280)
    assert L2_data.shape == (1280, 1280)

    assert L1b_nc_dqf.shape == (1280, 1280)
    assert L1b_fits_dqf.shape == (1280, 1280)
