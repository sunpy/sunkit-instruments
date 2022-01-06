import shutil
import os.path
import numpy
import pytest
import requests
import astropy
import sunpy
from sunkit_instruments import suvi
from sunkit_instruments.data.test import rootdir


@pytest.mark.remote_data
def test_suvi_download_from_NOAA(tmp_path):
    suvi.download_data_from_NOAA('2021-02-03T12:30:00', outdir=str(tmp_path), verbose=True)
    outfile = tmp_path / 'OR_SUVI-L1b-Fe171_G16_s20210341229311_e20210341229321_c20210341229525.fits.gz'
    assert outfile.exists()


@pytest.mark.remote_data
def test_suvi_testdata_exists(tmp_path):
    dst1 = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")
    dst2 = os.path.join(rootdir, "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits")
    dst3 = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc")
    
    if not os.path.exists(dst1):
        try:
            f1 = requests.get("https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/"+
                              "suvi-l1b-fe171/2021/12/31/OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")
            open(os.path.join(str(tmp_path), "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz"), 'wb').write(f1.content)
        finally:
            f1.close()
        outfile1 = tmp_path / "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz"
        shutil.copyfile(str(outfile1), dst1)
        
    if not os.path.exists(dst2):
        try:
            f2 = requests.get("https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/"+
                              "suvi-l2-ci171/2021/12/31/dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits")
            open(os.path.join(str(tmp_path), "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits"), 'wb').write(f2.content)
        finally:
            f2.close()
        outfile2 = tmp_path / "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits"
        shutil.copyfile(str(outfile2), dst2)
        
    if not os.path.exists(dst3):
        try:
            f3 = requests.get("https://noaa-goes16.s3.amazonaws.com/SUVI-L1b-Fe171/2021/365/00/"+
                              "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc")
            open(os.path.join(str(tmp_path), "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc"), 'wb').write(f3.content)
        finally:
            f3.close()
        outfile3 = tmp_path / "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc"
        shutil.copyfile(str(outfile3), dst3)

    assert os.path.exists(dst1)
    assert os.path.exists(dst2)
    assert os.path.exists(dst3)


def test_suvi_fix_L1b_header():
    L1b_fits = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")
    header = suvi.fix_L1b_header(L1b_fits)

    assert isinstance(header, astropy.io.fits.header.Header)


def test_suvi_despiking():
    L1b_nc = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc")
    L1b_fits = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")

    L1b_nc_header, L1b_nc_data, L1b_nc_dqf = suvi.read_suvi(L1b_nc, return_DQF=True)
    L1b_fits_header, L1b_fits_data, L1b_fits_dqf = suvi.read_suvi(L1b_fits, return_DQF=True)

    despiked_L1b_nc_data = L1b_nc_data
    despiked_L1b_fits_data = L1b_fits_data

    despiked_L1b_nc_data = suvi.despike_L1b_image(L1b_nc)
    despiked_L1b_fits_data = suvi.despike_L1b_image(L1b_fits)

    assert numpy.array_equal(L1b_nc_data, despiked_L1b_nc_data) is False
    assert numpy.array_equal(L1b_fits_data, despiked_L1b_fits_data) is False


def test_files_to_map():
    L1b_nc = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc")
    L1b_fits = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")
    L2_composite = os.path.join(rootdir, "dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits")

    L1b_nc_map = suvi.files_to_map(L1b_nc)
    L1b_fits_map = suvi.files_to_map(L1b_fits)
    L2_map = suvi.files_to_map(L2_composite)

    assert isinstance(L1b_nc_map, sunpy.map.sources.suvi.SUVIMap)
    assert isinstance(L1b_fits_map, sunpy.map.sources.suvi.SUVIMap)
    assert isinstance(L2_map, sunpy.map.sources.suvi.SUVIMap)


def test_response_functions():
    L1b_nc = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc")
    L1b_fits = os.path.join(rootdir, "OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz")

    L1b_nc_response = suvi.get_response(L1b_nc)
    L1b_fits_response = suvi.get_response(L1b_fits)
    response_195 = suvi.get_response(195)

    assert L1b_nc_response['wavelength_channel'] == 171
    assert L1b_fits_response['wavelength_channel'] == 171
    assert response_195['wavelength_channel'] == 195
