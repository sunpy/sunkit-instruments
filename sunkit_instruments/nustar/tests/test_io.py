from unittest.mock import MagicMock, patch

import numpy as np

from sunkit_instruments.nustar.io import read_heasarc_arf, read_nustar_pha, read_heasarc_rmf


@patch("astropy.io.fits.open")
def test_read_nustar_pha(mock_open):
    channel = np.array([0, 1, 2])
    counts = np.array([10, 20, 30])
    livetime = 1.0
    hdul = MagicMock()
    hdul[1].data = {"channel": np.array([0, 1, 2]), "counts": np.array([10, 20, 30])}
    hdul[0].header = {"LIVETIME": 1.0}
    mock_open.return_value.__enter__.return_value = hdul
    data, header = read_nustar_pha("test.pha")
    for t, d in zip((channel, counts), ("channel", "counts")):
        assert np.all(t == data[d])
    assert np.all(livetime==header["LIVETIME"])


@patch("astropy.io.fits.open")
def test_read_heasarc_arf(mock_open):
    e_lo = np.array([0, 1, 2]) 
    e_hi = np.array([1, 2, 3]) 
    area = np.array([10, 20, 30]) 
    hdul = MagicMock()
    hdul.header = {"HDUCLAS2":"SPECRESP"}
    hdul.data = {"energ_lo": np.array([0, 1, 2]), 
                 "energ_hi": np.array([1, 2, 3]),
                 "specresp":np.array([10, 20, 30])}
    mock_open.return_value.__enter__.return_value = (hdul,)
    data = read_heasarc_arf("test.arf")
    for t, d in zip((e_lo, e_hi, area), ("energ_lo", "energ_hi","specresp")):
        assert np.all(t == data[d])


@patch("astropy.io.fits.open")
def test_read_heasarc_rmf(mock_open):
    chan = np.array([0, 1, 2, 3])
    e_min = np.array([0.5, 1, 1.5, 2])
    e_max = np.array([1, 1.5, 2, 2.5])
    e_lo = np.array([0, 1, 2]) 
    e_hi = np.array([1, 2, 3]) 
    n_grp = np.array([10, 20, 30]) 
    f_chan = np.array([4, 5, 6]) 
    n_chan = np.array([40, 50, 60]) 
    matrix = np.array([-4, 8, 92]) 
    hdul0 = MagicMock()
    hdul1 = MagicMock()
    hdul0.header = {"HDUCLAS2":"EBOUNDS"}
    hdul0.data = {"channel":np.array([0, 1, 2, 3]), 
                  "e_min":np.array([0.5, 1, 1.5, 2]),
                  "e_max":np.array([1, 1.5, 2, 2.5])}
    hdul1.header = {"HDUCLAS2":"RSP_MATRIX"}
    hdul1.data = {"energ_lo":np.array([0, 1, 2]), 
                  "energ_hi":np.array([1, 2, 3]),
                  "n_grp":np.array([10, 20, 30]),
                  "f_chan":np.array([4, 5, 6]),
                  "n_chan":np.array([40, 50, 60]),
                  "matrix":np.array([-4, 8, 92])}
    mock_open.return_value.__enter__.return_value = (hdul0, hdul1)
    cdata, pdata = read_heasarc_rmf("test.rmf")
    for t, d in zip((chan, e_min, e_max), ("channel", "e_min","e_max")):
        assert np.all(t == cdata[d])
    for t, d in zip((e_lo, e_hi, n_grp, f_chan, n_chan, matrix), ("energ_lo", "energ_hi","n_grp", "f_chan", "n_chan", "matrix")):
        assert np.all(t == pdata[d])
