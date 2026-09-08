from unittest.mock import MagicMock, patch

import numpy as np

from sunkit_instruments.nustar.io import read_arf, read_pha, read_rmf


@patch("astropy.io.fits.open")
def test_read_pha(mock_open):
    channel = np.array([0, 1, 2])
    counts = np.array([10, 20, 30])
    livetime = 1.0
    hdul = MagicMock()
    hdul[1].data = {"channel": np.array([0, 1, 2]), "counts": np.array([10, 20, 30])}
    hdul[0].header = {"LIVETIME": 1.0}
    mock_open.return_value.__enter__.return_value = hdul
    data, header = read_pha("test.pha")
    for t, d in zip((channel, counts), ("channel", "counts")):
        assert np.all(t == data[d])
    assert np.all(livetime==header["LIVETIME"])


@patch("astropy.io.fits.open")
def test_read_arf(mock_open):
    e_lo = np.array([0, 1, 2]) 
    e_hi = np.array([1, 2, 3]) 
    area = np.array([10, 20, 30]) 
    hdul = MagicMock()
    hdul[1].data = {"energ_lo": np.array([0, 1, 2]), 
                    "energ_hi": np.array([1, 2, 3]),
                    "specresp":np.array([10, 20, 30])}
    mock_open.return_value.__enter__.return_value = hdul
    data = read_arf("test.arf")
    for t, d in zip((e_lo, e_hi, area), ("energ_lo", "energ_hi","specresp")):
        assert np.all(t == data[d])


@patch("astropy.io.fits.open")
def test_read_rmf(mock_open):
    e_lo = np.array([0, 1, 2]) 
    e_hi = np.array([1, 2, 3]) 
    n_grp = np.array([10, 20, 30]) 
    f_chan = np.array([4, 5, 6]) 
    n_chan = np.array([40, 50, 60]) 
    matrix = np.array([-4, 8, 92]) 
    hdul = MagicMock()
    hdul[2].data = {"energ_lo": np.array([0, 1, 2]), 
                    "energ_hi": np.array([1, 2, 3]),
                    "n_grp":np.array([10, 20, 30]),
                    "f_chan":np.array([4, 5, 6]),
                    "n_chan":np.array([40, 50, 60]),
                    "matrix":np.array([-4, 8, 92])}
    mock_open.return_value.__enter__.return_value = hdul
    data = read_rmf("test.rmf")
    for t, d in zip((e_lo, e_hi, n_grp, f_chan, n_chan, matrix), ("energ_lo", "energ_hi","n_grp", "f_chan", "n_chan", "matrix")):
        assert np.all(t == data[d])
