import warnings

import numpy as np

import astropy.units as u

from ndcube import NDMeta
from sunkit_spex.spectrum.spectrum import SpectralAxis, Spectrum
from sunkit_spex.spectrum.uncertainty import PoissonUncertainty

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
    """Store a bit of set-up values to reuse."""
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
    expected_rmf = np.array([[10.,  0., 20.,  0.],
                             [30.,  0.,  0.,  0.],
                             [ 0., 40.,  0.,  0.],
                             [ 0.,  0., 15., 16.]]) << (u.ct/u.ph)
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
            "data":data,
            "expected_rmf":expected_rmf,
            }

def _NustarSpectrum_setup0():
    """Return a `NustarSpectrum` object with the file functions re-defined."""
    setup = _NustarSpectrum_inputs_setup0()
    def obs_func(*args):
        return (setup["chan"], setup["counts"], setup["lvt"])
    def eff_func(*args): 
        return (setup["arf_elo"], setup["arf_ehi"], setup["arf_resp"])
    def res_func(*args): 
        return ((setup["rmf_chan"], setup["rmf_emin"], setup["rmf_emax"]), (setup["rmf_elo"], setup["rmf_ehi"], setup["ngrp"], setup["fchan"], setup["nchan"], setup["data"]))

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
    assert np.all(nu_spec._redistribution_matrix_ouput_channel_number==setup["rmf_chan"])
    assert np.all(nu_spec._redistribution_matrix_input_axis_edges==rmf_input_edges)
    assert np.all(nu_spec._redistribution_matrix_output_axis_edges==rmf_output_edges)
    for (k, v) in rmf_aux_info.items():
        assert np.all(nu_spec._redistribution_matrix_aux_info[k]==v)
    assert np.all(nu_spec._redistribution_matrix==setup["expected_rmf"])
    # SRM construction and assignment
    expected_srm = (setup["arf_resp"][:, None] * setup["expected_rmf"])
    assert np.all(nu_spec._spectral_response_matrix_input_axis_edges==rmf_input_edges)
    assert np.all(nu_spec._spectral_response_matrix_output_axis_edges==rmf_output_edges)
    assert np.all(nu_spec._spectral_response_matrix==expected_srm)

def test_NustarSpectrum_get_functions():
    """Test all `NustarSpectrum` functions that get data."""
    setup = _NustarSpectrum_inputs_setup0()
    nu_spec = _NustarSpectrum_setup0()

    rmf_input_edges = np.hstack((setup["rmf_elo"][:,None], setup["rmf_ehi"][:,None]))
    rmf_output_edges = np.hstack((setup["rmf_emin"][:,None], setup["rmf_emax"][:,None]))
    arf_edges = np.hstack((setup["arf_elo"][:,None], setup["arf_ehi"][:,None]))

    # check PHA information returned is correct
    pha_dict = nu_spec.get_pha_info()
    # check ARF information returned is correct
    arf_dict = nu_spec.get_arf_info()
    # check RMF information returned is correct
    rmf_dict = nu_spec.get_rmf_info()
    rmf_dict_aux = nu_spec.get_rmf_info(include_auxilliray_info=True)
    # check SRM information returned is correct
    srm_dict = nu_spec.get_srm_info()

    # make expected dictionaries
    expected_pha_dict = {"spectrum_counts_axis":rmf_output_edges,
                         "spectrum_counts":setup["counts"],
                         "effective_exposure":setup["lvt"]}
    expected_arf_dict = {"effective_area_axis":arf_edges,
                         "effective_area":setup["arf_resp"]}
    expected_rmf_dict = {"rmf_input_axis":rmf_input_edges,
                         "rmf_output_axis":rmf_output_edges,
                         "rmf":setup["expected_rmf"]}
    expected_rmf_dict_aux = expected_rmf_dict | {"e_lo_rmf":setup["rmf_elo"],
                                                 "e_hi_rmf":setup["rmf_ehi"],
                                                 "ngrp":setup["ngrp"],
                                                 "fchan":setup["fchan"],
                                                 "nchan":setup["nchan"],
                                                 "matrix":setup["data"]}
    expected_srm_dict = {"srm_input_axis":rmf_input_edges,
                         "srm_output_axis":rmf_output_edges,
                         "srm":(setup["arf_resp"][:, None] * setup["expected_rmf"])}

    for vd, ve in zip(pha_dict.values(), expected_pha_dict.values()):
        assert np.all(vd==ve)
    for vd, ve in zip(arf_dict.values(), expected_arf_dict.values()):
        assert np.all(vd==ve)
    for vd, ve in zip(rmf_dict.values(), expected_rmf_dict.values()):
        assert np.all(vd==ve)
    for vd, ve in zip(rmf_dict_aux.values(), expected_rmf_dict_aux.values()):
        assert np.all(vd==ve)
    for vd, ve in zip(srm_dict.values(), expected_srm_dict.values()):
        assert np.all(vd==ve)

def test_NustarSpectrum_spectrum_object():
    """Test `~NustarSpectrum.spectrum_object`."""
    setup = _NustarSpectrum_inputs_setup0()
    nu_spec = _NustarSpectrum_setup0()

    spec_obj_get = nu_spec.get_spec_obj()
    spec_obj_att = nu_spec.spectrum_object

    counts = setup["counts"]
    counts_uncertainity_pu = PoissonUncertainty(np.sqrt(counts.value)<<counts.unit)
    counts_spectral_axis = SpectralAxis(np.append(setup["rmf_emin"], setup["rmf_emax"][-1]), bin_specification="edges")
    meta = NDMeta()
    meta.add("exposure_time", setup["lvt"])
    meta.add("srm", setup["arf_resp"][:, None] * setup["expected_rmf"])
    meta.add("ph_axis", np.hstack((setup["rmf_elo"][:,None], setup["rmf_ehi"][:,None])))
    meta.add("distance", 1.0<<u.AU)
    spec_obj_made = Spectrum(
        data=counts, uncertainty=counts_uncertainity_pu, spectral_axis=counts_spectral_axis, meta=meta
    )

    assert np.all(spec_obj_get.data==spec_obj_att.data)
    assert np.all(spec_obj_get.data==spec_obj_made.data)
    assert np.all(spec_obj_get.uncertainty.array==spec_obj_att.uncertainty.array)
    assert np.all(spec_obj_get.uncertainty.array==spec_obj_made.uncertainty.array)
    assert np.all(spec_obj_get.spectral_axis==spec_obj_att.spectral_axis)
    assert np.all(spec_obj_get.spectral_axis==spec_obj_made.spectral_axis)
    for (k,v) in spec_obj_made.meta.items():
        assert np.all(v==spec_obj_get.meta[k])
        assert np.all(v==spec_obj_att.meta[k])
    assert np.all(list(spec_obj_get.meta.keys())==list(spec_obj_att.meta.keys()))
    assert np.all(list(spec_obj_get.meta.keys())==list(spec_obj_made.meta.keys()))
    
def test_NustarSpectrum_rebin_functions():
    """Test all `NustarSpectrum` functions that rebin data."""
    setup = _NustarSpectrum_inputs_setup0()
    nu_spec = _NustarSpectrum_setup0()
    # let's just turn everything into one bin
    _new_input_axis_edges = np.array([[setup["arf_elo"][0].value, setup["arf_ehi"][-1].value]]) << u.keV
    _new_output_axis_edges = np.array([[setup["rmf_emin"][0].value, setup["rmf_emax"][-1].value]]) << u.keV

    # check the individual rebin functions that only return the new arrays
    with warnings.catch_warnings(action="ignore"):
        # PHA rebinning
        new_axis_edges_pha, new_pha = nu_spec.rebin_pha_info(new_axis_edges=_new_output_axis_edges)
        # ARF rebinning
        new_axis_edges_arf, new_arf = nu_spec.rebin_arf_info(new_axis_edges=_new_input_axis_edges)
        # RMF rebinning
        new_input_axis_edges_rmf, new_output_axis_edges_rmf, new_rmf = nu_spec.rebin_rmf_info(new_input_axis_edges=_new_input_axis_edges, new_output_axis_edges=_new_output_axis_edges)
        # SRM rebinning, this rebins ARF and RMF separately then recombines
        new_input_axis_edges_srm, new_output_axis_edges_srm, new_arf_srm, new_rmf_srm, new_srm = nu_spec.rebin_srm_info(new_input_axis_edges=_new_input_axis_edges, new_output_axis_edges=_new_output_axis_edges)

    assert np.all(new_axis_edges_pha==_new_output_axis_edges)
    assert np.all(new_pha==np.sum(setup["counts"]))
    assert np.all(new_axis_edges_arf==_new_input_axis_edges)
    assert np.all(new_arf==np.mean(setup["arf_resp"]))
    assert np.all(new_input_axis_edges_rmf==_new_input_axis_edges)
    assert np.all(new_output_axis_edges_rmf==_new_output_axis_edges)
    assert np.all(new_rmf==np.sum(np.mean(setup["expected_rmf"], axis=0))) # mean over rows, sum columns

    expected_srm = (new_arf_srm[:, None] * new_rmf_srm)
    assert np.all(new_input_axis_edges_srm==_new_input_axis_edges)
    assert np.all(new_output_axis_edges_srm==_new_output_axis_edges)
    assert np.all(new_arf==new_arf_srm)
    assert np.all(new_rmf==new_rmf_srm)
    assert np.all(new_srm==expected_srm)

    # rebin all info; this sets the class attributes and does not return anything
    nu_spec.rebin_info(new_input_axis_edges=_new_input_axis_edges, new_output_axis_edges=_new_output_axis_edges)
    pha_dict = nu_spec.get_pha_info()
    arf_dict = nu_spec.get_arf_info()
    rmf_dict = nu_spec.get_rmf_info()
    srm_dict = nu_spec.get_srm_info()
    for (v,e) in zip(pha_dict.values(), [new_axis_edges_pha, new_pha, setup["lvt"]]):
        assert np.all(v==e)
    for (v,e) in zip(arf_dict.values(), [new_axis_edges_arf, new_arf]):
        assert np.all(v==e)
    for (v,e) in zip(rmf_dict.values(), [new_input_axis_edges_rmf, new_output_axis_edges_rmf, new_rmf]):
        assert np.all(v==e)
    for (v,e) in zip(srm_dict.values(), [new_input_axis_edges_srm, new_output_axis_edges_srm, new_srm]):
        assert np.all(v==e)

def test_NustarSpectrum__define_standard_units():
    """Test `~NustarSpectrum._define_standard_units` function."""
    nu_spec = _NustarSpectrum_setup0()
    standard_units = {"channel_number":(u.dimensionless_unscaled),
                      "energy":(u.keV),
                      "ct_spec":(u.ct),
                      "eff_exp/lvt":(u.s),
                      "eff_area":(u.cm**2),
                      "rdm":(u.ct * u.ph**-1),
                      "srm":(u.ct * u.ph**-1 * u.cm**2)}
    for ((k,v),(ek,ev)) in zip(nu_spec.get_standard_units().items(), standard_units.items()):
        assert k==ek
        assert v==ev
