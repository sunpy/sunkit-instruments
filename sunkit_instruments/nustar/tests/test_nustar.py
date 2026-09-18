import warnings

import numpy as np

import astropy.units as u

from sunkit_instruments.nustar.nustar import NustarSpectrum


def _init_NustarSpectrum(obs_func, eff_func, res_func):
    """For testing, setup and redefine some methods that rely on files."""

    class NonFileNustarSpectrum(NustarSpectrum):
        def __init__(self):
            super().__init__("test.pha",
                             arf_file="test.arf",
                             rmf_file="test.rmf")
        def _get_observable_info(self, *args):
            return obs_func(*args)
        def _get_effective_area_info(self, *args):
            return eff_func(*args)
        def _get_response_info(self, *args):
            return res_func(*args)

    with warnings.catch_warnings(action="ignore"):
        return NonFileNustarSpectrum()

def _NustarSpectrum_inputs_setup0():
    """Store a but of set-up values to reuse."""
    chan = np.array([0, 1, 2, 3]) << u.dimensionless_unscaled
    counts = np.array([8, 4, 9, 2.3]) << u.ct
    lvt = 11 << u.second

    arf_elo = np.array([1.5, 2.5, 3.5, 4.5]) << u.keV
    arf_ehi = np.array([2.5, 3.5, 4.5, 5.5]) << u.keV
    arf_resp = np.array([14, 12, 7, 3]) << u.cm**2

    rmf_chan = np.array([0, 1, 2, 3]) << u.dimensionless_unscaled
    rmf_emin = np.array([1.5, 2.5, 3.5, 4.5]) << u.keV
    rmf_emax = np.array([2.5, 3.5, 4.5, 5.5]) << u.keV
    rmf_elo = np.array([1.5, 2.5, 3.5, 4.5]) << u.keV
    rmf_ehi = np.array([2.5, 3.5, 4.5, 5.5]) << u.keV
    ngrp = np.array([2, 1, 1, 2]) << u.dimensionless_unscaled
    fchan = [[0, 2], [0], [1], [2]]
    nchan = [[1, 1], [1], [1], [2]]
    data = [[10, 20], [30], [40], [15, 16]]
    return {"chan":chan,
            "counts":counts,
            "lvt":lvt,
            "arf_elo":arf_elo,
            "arf_ehi":arf_ehi,
            "arf_resp":arf_resp,
            "rmf_chan":rmf_chan,
            "rmf_emin":rmf_emin,
            "rmf_emax":rmf_emax,
            "rmf_elo":rmf_elo,
            "rmf_ehi":rmf_ehi,
            "ngrp":ngrp,
            "fchan":fchan,
            "nchan":nchan,
            "data":data}

def _NustarSpectrum_setup0():
    """Return a `NustarSpectrum` object with the file functions re-defined."""
    setup = _NustarSpectrum_inputs_setup0()
    obs_func = lambda *args: (setup["chan"], setup["counts"], setup["lvt"])
    eff_func = lambda *args: (setup["arf_elo"], setup["arf_ehi"], setup["arf_resp"])
    res_func = lambda *args: ((setup["rmf_chan"], setup["rmf_emin"], setup["rmf_emax"]), (setup["rmf_elo"], setup["rmf_ehi"], setup["ngrp"], setup["fchan"], setup["nchan"], setup["data"]))

    return _init_NustarSpectrum(obs_func, eff_func, res_func)

def test_NustarSpectrum_assignment():
    """Test the first the data assignment in the `NustarSpectrum` class."""
    setup = _NustarSpectrum_inputs_setup0()
    nu_spec = _NustarSpectrum_setup0()

    # PHA assignment
    rmf_output_edges = np.hstack((setup["rmf_emin"][:,None], setup["rmf_emax"][:,None]))
    assert np.all(nu_spec._spectrum_channel_number==setup["chan"])
    assert np.all(nu_spec._spectrum_counts==setup["counts"])
    assert np.all(nu_spec._effective_exposure==setup["lvt"])
    assert np.all(nu_spec._spectrum_axis_edges==rmf_output_edges)
    # ARF assignment
    arf_edges = np.hstack((setup["arf_elo"][:,None], setup["arf_ehi"][:,None]))
    assert np.all(nu_spec._effective_area==setup["arf_resp"])
    assert np.all(nu_spec._effective_area_axis_edges==arf_edges)
    # RMF construction and assignment
    rmf_input_edges = np.hstack((setup["rmf_elo"][:,None], setup["rmf_ehi"][:,None]))
    rmf_aux_info = {"e_lo_rmf":setup["rmf_elo"],
                    "e_hi_rmf":setup["rmf_ehi"],
                    "ngrp":setup["ngrp"],
                    "fchan":setup["fchan"],
                    "nchan":setup["nchan"],
                    "matrix":setup["data"]}
    expected_rmf = np.array([[10.,  0., 20.,  0.],
                             [30.,  0.,  0.,  0.],
                             [ 0., 40.,  0.,  0.],
                             [ 0.,  0., 15., 16.]]) << (u.ct/u.ph)
    assert np.all(nu_spec._redistribution_matrix_ouput_channel_number==setup["rmf_chan"])
    assert np.all(nu_spec._redistribution_matrix_input_axis_edges==rmf_input_edges)
    assert np.all(nu_spec._redistribution_matrix_output_axis_edges==rmf_output_edges)
    for (k, v) in rmf_aux_info.items():
        assert np.all(nu_spec._redistribution_matrix_aux_info[k]==v)
    assert np.all(nu_spec._redistribution_matrix==expected_rmf)
    # SRM construction and assignment
    expected_srm = (setup["arf_resp"][:, None] * expected_rmf)
    assert np.all(nu_spec._spectral_response_matrix_input_axis_edges==rmf_input_edges)
    assert np.all(nu_spec._spectral_response_matrix_output_axis_edges==rmf_output_edges)
    assert np.all(nu_spec._spectral_response_matrix==expected_srm)

def test_NustarSpectrum_spectrum_object():
    """Test `~NustarSpectrum.spectrum_object`."""

def test_NustarSpectrum_get_functions():
    """Test all `NustarSpectrum` functions that get data."""

def test_NustarSpectrum_rebin_functions():
    """Test all `NustarSpectrum` functions that rebin data."""
