from tempfile import TemporaryDirectory

import pytest
from parfive import Downloader


@pytest.fixture(scope="session")
@pytest.mark.remote_data
def L1B_FITS():
    downloader = Downloader()
    url = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/suvi-l1b-fe171/2021/12/31/OR_SUVI-L1b-Fe171_G16_s20213650006108_e20213650006118_c20213650006321.fits.gz"
    with TemporaryDirectory() as d:
        downloader.enqueue_file(url, d)
        files = downloader.download()
        yield files[0]


@pytest.fixture(scope="session")
@pytest.mark.remote_data
def L1B_NC():
    downloader = Downloader()
    url = "https://noaa-goes16.s3.amazonaws.com/SUVI-L1b-Fe171/2021/365/00/OR_SUVI-L1b-Fe171_G16_s20213650022109_e20213650022119_c20213650022323.nc"
    with TemporaryDirectory() as d:
        downloader.enqueue_file(url, d)
        files = downloader.download()
        yield files[0]


@pytest.fixture(scope="session")
@pytest.mark.remote_data
def L2_COMPOSITE():
    downloader = Downloader()
    url = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/suvi-l2-ci171/2021/12/31/dr_suvi-l2-ci171_g16_s20211231T000800Z_e20211231T001200Z_v1-0-1.fits"
    with TemporaryDirectory() as d:
        downloader.enqueue_file(url, d)
        files = downloader.download()
        yield files[0]
