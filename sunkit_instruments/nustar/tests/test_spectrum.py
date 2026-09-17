from unittest.mock import MagicMock

import astropy.units as u
import numpy as np

from sunkit_instruments.nustar.spectrum import (
    get_observable_info, 
    get_effective_area_info, 
    col2arr, 
    vrmf2arr, 
    get_response_info, 
    make_srm,
    )


def test_get_observable_info():
    channel = np.array([0, 1, 2])
    counts = np.array([10, 20, 30])
    livetime = 1.0
    hdul = MagicMock()
    hdul[1].data = {"channel": np.array([0, 1, 2]), "counts": np.array([10, 20, 30])}
    hdul[0].header = {"LIVETIME": 1.0}
    ch, ct, lvt = get_observable_info(hdul[1].data, hdul[0].header)
    assert np.all((channel<<u.dimensionless_unscaled)==ch)
    assert np.all((counts<<u.ct)==ct)
    assert np.all((livetime<<u.s)==lvt)

def test_get_effective_area_info():
    e_lo = np.array([0, 1, 2]) 
    e_hi = np.array([1, 2, 3]) 
    area = np.array([10, 20, 30]) 
    hdul = MagicMock()
    hdul[1].data = {"energ_lo": np.array([0, 1, 2]), 
                    "energ_hi": np.array([1, 2, 3]),
                    "specresp":np.array([10, 20, 30])}
    el, eh, a = get_effective_area_info(hdul[1].data)
    assert np.all((e_lo<<u.keV)==el)
    assert np.all((e_hi<<u.keV)==eh)
    assert np.all((area<<u.cm**2)==a)

def test_get_response_info():
    chan = np.array([0, 1, 2, 3])
    e_min = np.array([0.5, 1, 1.5, 2])
    e_max = np.array([1, 1.5, 2, 2.5])
    e_lo = np.array([0, 1, 2]) 
    e_hi = np.array([1, 2, 3]) 
    n_grp = np.array([10, 20, 30]) 
    f_chan = np.array([4, 5, 6]) 
    n_chan = np.array([40, 50, 60]) 
    matrix = np.array([-4, 8, 92]) 
    cdata = {"channel":np.array([0, 1, 2, 3]), 
             "e_min":np.array([0.5, 1, 1.5, 2]),
             "e_max":np.array([1, 1.5, 2, 2.5])}
    pdata = {"energ_lo": np.array([0, 1, 2]), 
             "energ_hi": np.array([1, 2, 3]),
             "n_grp":np.array([10, 20, 30]),
             "f_chan":np.array([4, 5, 6]),
             "n_chan":np.array([40, 50, 60]),
             "matrix":np.array([-4, 8, 92])}
    (c, emi, ema), (el, eh, ng, fc, nc, m) = get_response_info(cdata, pdata)
    assert np.all((chan<<u.dimensionless_unscaled)==c)
    assert np.all((e_min<<u.keV)==emi)
    assert np.all((e_max<<u.keV)==ema)
    assert np.all((e_lo<<u.keV)==el)
    assert np.all((e_hi<<u.keV)==eh)
    assert np.all((n_grp<<u.dimensionless_unscaled)==ng)
    assert np.all((f_chan<<u.dimensionless_unscaled)==fc)
    assert np.all((n_chan<<u.dimensionless_unscaled)==nc)
    assert np.all(matrix==m)

def test_col2arr():
    row_data = [[0],
                [0],
                [0, 22],
                [0, 12, 40]]
    expected = np.array([[0, 0, 0],
                         [0, 0 ,0],
                         [0, 22, 0],
                         [0, 12, 40]])
    output_array = col2arr(row_data)
    assert np.all(expected==output_array)

def test_vrmf2arr():
    data = [[10, 20], [30], [40]]
    ngrp = np.array([2, 1, 1])
    fchan = np.array([[0, 2], [0, 0], [1, 0]])
    nchan = np.array([[1, 1], [1, 0], [1, 0]])
    output_array = vrmf2arr(data=data, n_grp_list=ngrp, f_chan_array=fchan, n_chan_array=nchan)
    expected = np.array([[10.,  0., 20.],
                         [30.,  0.,  0.],
                         [ 0., 40.,  0.]]) << (u.ct / u.ph)
    assert np.all(expected==output_array)

def test_make_srm():
    rmf_matrix = np.array([[10.,  0., 20.],
                           [30.,  0.,  0.],
                           [ 0., 40.,  0.],
                           [ 0., 50.,  60.]]) << (u.ct / u.ph)
    arf_array = np.array([2, 0.5, 0.1, 0.01]) << u.cm**2
    output_array = make_srm(rmf_matrix, arf_array)
    expected = np.array([[20.,  0., 40.],
                         [15.,  0.,  0.],
                         [ 0.,  4.,  0.],
                         [ 0., 0.5,  0.6]]) << (u.ct * u.cm**2 / u.ph)
    assert np.all(expected==output_array)
