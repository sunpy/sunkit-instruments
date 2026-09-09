import os
import warnings

import astropy.units as u
import numpy as np

from sunkit_instruments.nustar.io import (
    read_arf,
    read_pha,
    read_rmf,
)
from sunkit_instruments.nustar.spectrum import (
    get_observable_info,
    get_effective_area_info,
    get_response_info,
    col2arr,
    vrmf2arr,
    make_srm,
)
from sunkit_instruments.nustar.utils import (
    rebin_rmf,
    regroup_any_array,
)

class NustarSpectrum:
    """
    Loader specifically for NuSTAR spectral files.

    Provides a nice-to-have class for loading in NuSTAR `.pha`, `.arf`, 
    and `.rmf` files already produced using the NuSTARDAS software. This
    data can then be inspected using the class and other Python tools.

    The method `get_spec_obj` will then be able to return a spectrum 
    data container to be used in spectral fitting software tools like
    Python's `sunkit_spec` [1].

    [1] https://sunkit-spex.readthedocs.io/en/latest/

    Parameters
    ----------
    pha_file : `str`
            The PHA file of the spectrum to be loaded.

    arf_file, rmf_file : `str`, `str`
            The ARF and RMF files associated with the PHA file(s). 
            Default: None, None
    """

    def __init__(self, pha_file, arf_file=None, rmf_file=None):
        self._standard_units = {"channel_number":(u.dimensionless_unscaled),
                                "energy":(u.keV),
                                "ct_spec":(u.ct),
                                "eff_exp/lvt":(u.s),
                                "eff_area":(u.cm**2),
                                "rdm":(u.ct * u.ph**-1),
                                "srm":(u.ct * u.ph**-1 * u.cm**2)}
        self._construction_string = (
            f"NustarSpectrum(pha_file={pha_file},arf_file={arf_file},rmf_file={rmf_file})"
        )
        # RMF contains channel->count energy information so do this first if we can
        self.construct_rmf(f_rmf=rmf_file)
        self.construct_spectrum(pha_file)
        self.construct_arf(f_arf=arf_file)
        # the SRM construction needs the RMF and ARF to be done
        self.construct_srm()

    def get_spec_obj(self):
        """..."""
        pass

    def construct_spectrum(self, f_pha):
        self._spectrum_channel_number, self._spectrum_counts, self._effective_exposure = get_observable_info(*read_pha(f_pha))

        self.standard_unit_check("Count spectrum channel numbers", 
                                 self._spectrum_channel_number.unit, 
                                 "channel_number")
        self.standard_unit_check("Count spectrum array", 
                                 self._spectrum_counts.unit, 
                                 "ct_spec")
        self.standard_unit_check("Effective exposure/livetime value", 
                                 self._effective_exposure.unit, 
                                 "eff_exp/lvt")

        # if an RMF exists then we can map the channel number to energy, else just assume standard for convenience
        if self._has_rmf:
            indices = [np.nonzero(self._redistribution_matrix_ouput_channel_number.value==chan)[0] for chan in self._spectrum_channel_number]
            self._spectrum_axis_edges = self._redistribution_matrix_output_axis_edges[indices,:].squeeze()
            if not np.allclose(self._spectrum_axis_edges, self._redistribution_matrix_output_axis_edges):
                warnings.warn("Spectrum and RMF count bin edge information is different.")
        else:
            warnings.warn(f"No RMF information exists so defaulting to standard NuSTAR energy binning for the spectrum.")
            self.set_standard_spectrum_axis()

    def set_standard_spectrum_axis(self):
        _standard_energy_start = 1.6 << u.keV
        _standard_num_of_channels = 4096
        _standard_energy_binning = 0.04 << u.keV
        _standard_max_energy = _standard_num_of_channels*_standard_energy_binning + _standard_energy_start
        e_lo = np.arange(_standard_energy_start.value, _standard_max_energy.value, _standard_energy_binning.value) << _standard_energy_start.unit
        e_hi = e_lo + _standard_energy_binning
        self._spectrum_axis_edges = np.hstack((e_lo[:,None], e_hi[:,None]))
        
    def construct_arf(self, f_arf=None):
        if (f_arf is None) or (not os.path.isfile(f_arf)):
            warnings.warn(f"File `{f_arf}` is not found or has not been given.")
            self._effective_area_axis_edges = None
            self._effective_area = None
            self._has_arf = False
            return

        e_lo_arf, e_hi_arf, self._effective_area = get_effective_area_info(read_arf(f_arf))
        self._effective_area_axis_edges = np.hstack((e_lo_arf[:,None], e_hi_arf[:,None]))
        self._has_arf = True

        self.standard_unit_check("Effective area array", 
                                 self._effective_area.unit, 
                                 "eff_area")
        self.standard_unit_check("Effective area axis", 
                                 self._effective_area_axis_edges.unit, 
                                 "energy")

    def construct_rmf(self, f_rmf=None):
        if (f_rmf is None) or (not os.path.isfile(f_rmf)):
            warnings.warn(f"File `{f_rmf}` is not found or has not been given.")
            self._redistribution_matrix_aux_info = None
            self._redistribution_matrix_input_axis_edges = None
            self._redistribution_matrix_ouput_channel_number = None
            self._redistribution_matrix_output_axis_edges = None
            self._redistribution_matrix = None
            self._has_rmf = False
            return
        
        (chan, e_min, e_max), (e_lo_rmf, e_hi_rmf, ngrp, fchan, nchan, matrix, redist_m) = self.load_rmf(f_rmf)

        self._redistribution_matrix_aux_info = {"e_lo_rmf":e_lo_rmf, 
                                                "e_hi_rmf":e_hi_rmf, 
                                                "ngrp":ngrp, 
                                                "fchan":fchan, 
                                                "nchan":nchan, 
                                                "matrix":matrix}
        self._redistribution_matrix_input_axis_edges = np.hstack((e_lo_rmf[:,None], e_hi_rmf[:,None]))
        self._redistribution_matrix_ouput_channel_number = chan
        self._redistribution_matrix_output_axis_edges = np.hstack((e_min[:,None],  e_max[:,None]))
        self._redistribution_matrix = redist_m
        self._has_rmf = True

        self.standard_unit_check("RMF count channel numbers", 
                                 self._redistribution_matrix_ouput_channel_number.unit, 
                                 "channel_number")
        self.standard_unit_check("Redistribution matrix", 
                                 self._redistribution_matrix.unit, 
                                 "rdm")
        self.standard_unit_check("Redistribution matrix input axis", 
                                 self._redistribution_matrix_input_axis_edges.unit, 
                                 "energy")
        self.standard_unit_check("Redistribution matrix output axis", 
                                 self._redistribution_matrix_output_axis_edges.unit, 
                                 "energy")

    def construct_srm(self):
        """..."""
        self._has_srm = False
        if not self._has_arf:
            warnings.warn("Cannot construct SRM as ARF information is missing.")
            return
        if not self._has_rmf:
            warnings.warn("Cannot construct SRM as RMF information is missing.")
            return
        if not np.allclose(self._effective_area_axis_edges, self._redistribution_matrix_input_axis_edges):
            warnings.warn("RMF and ARF information are formatted for different axes.")
            return

        self._spectral_response_matrix = make_srm(rmf_matrix=self._redistribution_matrix, 
                                                  arf_array=self._effective_area)
        self._has_srm = True
        
        self.standard_unit_check("Spectral response matrix", 
                                 self._spectral_response_matrix.unit, 
                                 "srm")
            
    def get_srm(self):
        """..."""
        if not self._has_srm:
            warnings.warn("Missing SRM information.")
            return
        return self._spectral_response_matrix

    def get_count_energy_bins(self):
        """..."""
        if not self._has_rmf:
            warnings.warn("Missing count energy bins as no RMF information is available. Returning default spectrum axes instead.")
            return self._spectrum_axis_edges

        if not np.allclose(self._spectrum_axis_edges, self._redistribution_matrix_output_axis_edges):
            warnings.warn("Spectrum and RMF count bin edge information is different. Returning dictionary of both.")
            return {"spec-count-bin-edges":self._spectrum_axis_edges,
                    "rmf-count-bin-edges":self._redistribution_matrix_input_axis_edges}
        
        return self._redistribution_matrix_output_axis_edges

    def get_photon_energy_bins(self):
        """..."""
        if (not self._has_arf) and (self._has_rmf):
            warnings.warn("Missing photon energy bin information as ARF is not available. Only returning information from RMF.")
            return self._redistribution_matrix_input_axis_edges
        if (self._has_arf) and (not self._has_rmf):
            warnings.warn("Missing photon energy bin information as RMF is not available. Only returning information from ARF.")
            return self._effective_area_axis_edges
        if (not self._has_arf) and (not self._has_rmf):
            warnings.warn("Missing photon energy bin information as ARF and RMF is not available")
            return 
        
        if not np.allclose(self._effective_area_axis_edges, self._redistribution_matrix_input_axis_edges):
            warnings.warn("RMF and ARF photon bin edge information is different. Returning dictionary of both.")
            return {"arf-photon-bin-edges":self._effective_area_axis_edges,
                    "rmf-photon-bin-edges":self._redistribution_matrix_input_axis_edges}
        
        # if everything is fine then just return one, simple, nice array
        return self._redistribution_matrix_input_axis_edges

    def standard_unit_check(self, name, unit, key):
        if unit!=self._standard_units[key]:
            warnings.warn(f"{name} units are not standard {self._standard_units[key]}, but in {unit}.")

    def load_rmf(self, rmf_file):
        """Extracts all information, mainly the redistribution matrix ([counts/photon]) from a given RMF file.

        Parameters
        ----------
        rmf_file : string
                The file path and name of the RMF file.

        Returns
        -------
        The lower/higher photon bin edges (e_lo_rmf, e_hi_rmf), the number of counts channels activated by each photon channel (ngrp),
        starting indices of the count channel groups (fchan), number counts channels from each starting index (nchan), the corresponding
        counts/photon value for each count and photon entry (matrix), and the redistribution matrix (redist_m: with rows of photon channels,
        columns of counts channels, and in the units of counts/photon).
        """

        (chan, e_min, e_max), (e_lo_rmf, e_hi_rmf, ngrp, fchan, nchan, matrix) = get_response_info(*read_rmf(rmf_file))
        fchan_array = col2arr(fchan)
        nchan_array = col2arr(nchan)
        redist_m = vrmf2arr(
            data=matrix, n_grp_list=ngrp, f_chan_array=fchan_array, n_chan_array=nchan_array
        )  

        return (chan, e_min, e_max), (e_lo_rmf, e_hi_rmf, ngrp, fchan, nchan, matrix, redist_m)

    def _rebin_srm(self, axis="count"):
        """Rebins the photon and/or count channels of the spectral response matrix by rebinning the redistribution matrix and the effective area array.

        Parameters
        ----------
        axis : string
                Define what \'axis\' the binning should be applied to. E.g., \'photon\', \'count\', or \'photon_and_count\'.

        Returns
        -------
        The rebinned 2d spectral response matrix.
        """
        old_count_bins, new_count_bins, old_photon_bins, new_photon_bins = self._channel_bin_info(axis)

        old_rmf = self._loaded_spec_data["extras"]["rmf.redistribution_matrix"]
        old_eff_area = self._loaded_spec_data["extras"]["arf.effective_area"]

        # checked with ftrbnrmf
        new_rmf = rebin_rmf(
            matrix=old_rmf,
            old_count_bins=old_count_bins,
            new_count_bins=new_count_bins,
            old_photon_bins=old_photon_bins,
            new_photon_bins=new_photon_bins,
            axis=axis,
        )

        # average eff_area, checked with ftrbnarf
        new_eff_area = (
            regroup_any_array(data=old_eff_area, old_bins=old_photon_bins, new_bins=new_photon_bins, combine_by="mean")
            if (axis != "count")
            else old_eff_area
        )
        return make_srm(rmf_matrix=new_rmf, arf_array=new_eff_area)

    def __repr__(self):
        """String representation of `_loaded_spec_data`."""
        return self._construction_string
