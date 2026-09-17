import warnings

import astropy.units as u
import numpy as np

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
    """Store a but of set-up values to re-use."""
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

    assert np.all(nu_spec._spectrum_channel_number==setup["chan"])