import os
import warnings

import astropy.units as u
import numpy as np

from sunkit_instruments.nustar.io import (
    read_heasarc_arf,
    read_nustar_pha,
    read_heasarc_rmf,
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

__all__ = ["NustarSpectrum"]

class NustarSpectrum:
    """
    Loader specifically for NuSTAR spectral files.

    Provides a nice-to-have class for loading in NuSTAR `.pha`, `.arf`, 
    and `.rmf` files already produced using the NuSTARDAS software. This
    data can then be inspected using the class and other Python tools.

    The class is primarily built to be used with all 3 aforementioned
    files but will work with a subset of them. If a subset is given then 
    expect some useful warnings explaining the limitations.

    The method `get_spec_obj` will then be able to return a spectrum 
    data container to be used in spectral fitting software tools like
    Python's `sunkit_spec` [1].

    [1] https://sunkit-spex.readthedocs.io/en/latest/

    Parameters
    ----------
    pha_file : `str`
            The PHA file of the spectrum to be loaded.
            Default: None

    arf_file : `str`
            The ARF file, likely associated with the PHA file. 
            Default: None

    rmf_file : `str`
            The RMF file, likely associated with the PHA file. 
            Default: None, None
    """

    def __init__(self, pha_file=None, arf_file=None, rmf_file=None):
        self._standard_units = {"channel_number":(u.dimensionless_unscaled),
                                "energy":(u.keV),
                                "ct_spec":(u.ct),
                                "eff_exp/lvt":(u.s),
                                "eff_area":(u.cm**2),
                                "rdm":(u.ct * u.ph**-1),
                                "srm":(u.ct * u.ph**-1 * u.cm**2)}
        self._construction_string = (
            f"{str(self.__class__)}(pha_file={pha_file},arf_file={arf_file},rmf_file={rmf_file})"
        )
        # RMF contains channel->count energy information so do this first if we can
        self.construct_rmf(f_rmf=rmf_file)
        self.construct_pha(pha_file)
        self.construct_arf(f_arf=arf_file)
        # the SRM construction needs the RMF and ARF to be done
        self.construct_srm()

    def get_spec_obj(self):
        """..."""
        pass

    def construct_pha(self, f_pha=None):
        if (f_pha is None) or (not os.path.isfile(f_pha)):
            warnings.warn(f"File `{f_pha}` is not found or has not been given.")
            self._spectrum_channel_number = None
            self._spectrum_counts = None
            self._effective_exposure = None
            self._has_pha = False
            return

        self._spectrum_channel_number, self._spectrum_counts, self._effective_exposure = get_observable_info(*read_nustar_pha(f_pha))
        self._has_pha = True

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

        e_lo_arf, e_hi_arf, self._effective_area = get_effective_area_info(read_heasarc_arf(f_arf))
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
        if (not self._has_arf) and (not self._has_rmf):
            warnings.warn("Cannot construct SRM as ARF and RMF information is missing.")
            return
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
        self._spectral_response_matrix_input_axis_edges = self._redistribution_matrix_input_axis_edges
        self._spectral_response_matrix_output_axis_edges = self._redistribution_matrix_output_axis_edges
        self._has_srm = True
        
        self.standard_unit_check("Spectral response matrix", 
                                 self._spectral_response_matrix.unit, 
                                 "srm")
            
    def get_pha_info(self):
        """..."""
        if self._has_pha:
            return {"spectrum_counts_axis":self._spectrum_axis_edges,
                    "spectrum_counts":self._spectrum_counts}
        
        warnings.warn("Missing PHA information.")
            
    def set_pha_info(self, spectrum_counts=None, spectrum_counts_axis=None, pha_file=None):
        """..."""
        if pha_file is None:
            if spectrum_counts is not None:
                self._set_pha(spectrum_counts)
            if spectrum_counts_axis is not None:
                warnings.warn("""
                Editing this axis will affect the PHA's compatibility 
                with the RMF output axis.
                """
                )
                self._set_pha_axis(spectrum_counts_axis)
        else:
            self.construct_pha()

        self._post_pha_update_checks()

        if self._has_pha:
            warnings.warn("Overwritten PHA information.")
        self._has_pha = True

    def _post_pha_update_checks(self):
        """..."""
        # check the shapes match at least
        pha_shape = np.shape(self._spectrum_counts)
        if pha_shape[0]!=np.shape(self.spectrum_axis_edges)[0]:
            warnings.warn("Count axis length of PHA does not match PHA length.")

    def _set_pha(self, pha):
        """..."""
        self._spectrum_counts = pha 

        self.standard_unit_check("Count spectrum array", 
                                 self._spectrum_counts.unit, 
                                 "ct_spec")

    def _set_pha_axis(self, axis):
        """..."""
        self._axis_update("_spectrum_axis_edges", 
                          axis, 
                          "Spectrum energy axis")
            
    def get_arf_info(self):
        """..."""
        if self._has_arf:
            return {"effective_area_axis":self._effective_area_axis_edges,
                    "effective_area":self._effective_area}
        
        warnings.warn("Missing ARF information.")
            
    def set_arf_info(self, effective_area=None, effective_area_axis=None, arf_file=None):
        """..."""
        if arf_file is None:
            if effective_area is not None:
                self._set_arf(effective_area)
            if effective_area_axis is not None:
                warnings.warn("""
                Editing this axis will affect the ARF's compatibility 
                with the RMF input axis.
                """
                )
                self._set_arf_axis(effective_area_axis)
        else:
            self.construct_arf()

        self._post_arf_update_checks()

        if self._has_arf:
            warnings.warn("Overwritten ARF information.")
        self._has_arf = True

    def _post_arf_update_checks(self):
        """..."""
        # check the shapes match at least
        arf_shape = np.shape(self._effective_area)
        if arf_shape[0]!=np.shape(self._effective_area_axis_edges)[0]:
            warnings.warn("Photon axis length of ARF does not match ARF length.")

    def _set_arf(self, arf):
        """..."""
        self._effective_area = arf 

        self.standard_unit_check("Effective area array", 
                                 self._effective_area.unit, 
                                 "eff_area")

    def _set_arf_axis(self, axis):
        """..."""
        self._axis_update("_effective_area_axis_edges", 
                          axis, 
                          "Effective area axis")
            
    def get_rmf_info(self, include_auxilliray_info=False):
        """..."""
        if self._has_rmf:
            info = {"rmf_input_axis":self._redistribution_matrix_input_axis_edges,
                    "rmf_output_axis":self._redistribution_matrix_output_axis_edges,
                    "rmf":self._redistribution_matrix}

            if include_auxilliray_info:
                return info | self._redistribution_matrix_aux_info

            return info
        
        warnings.warn("Missing RMF information.")

    def set_rmf_info(self, rmf=None, output_axis_edges=None, input_axis_edges=None, rmf_file=None):
        """..."""

        if rmf_file is None:
            if rmf is not None:
                self._set_rmf(rmf)
            if output_axis_edges is not None:
                warnings.warn("""
                Editing the output axis will affect the RMF's compatibility 
                with the PHA, spectrum data.
                """
                )
                self._set_rmf_output_axis(output_axis_edges)
            if input_axis_edges is not None:
                warnings.warn("""
                Editing the input axis will affect the RMF's compatibility 
                with the ARF.
                """
                )
                self._set_rmf_input_axis(input_axis_edges)
        else:
            self.construct_rmf(f_rmf=rmf_file)

        self._post_rmf_update_checks()

        if self._has_rmf:
            warnings.warn("Overwritten RMF information.")
        self._has_rmf = True

    def _post_rmf_update_checks(self):
        """..."""
        # check the shapes match at least
        rmf_shape = np.shape(self._redistribution_matrix)
        if rmf_shape[0]!=np.shape(self._redistribution_matrix_input_axis_edges)[0]:
            warnings.warn("Input photon axis length of RMF does not match input axis values.")
        if rmf_shape[1]!=np.shape(self._redistribution_matrix_output_axis_edges)[0]:
            warnings.warn("Output count axis length of RMF does not match output axis values.")

    def _set_rmf(self, rmf):
        """..."""
        self._redistribution_matrix = rmf 

        self.standard_unit_check("Redistribution matrix", 
                                 self._redistribution_matrix.unit, 
                                 "rmf")

    def _set_rmf_input_axis(self, input_axis):
        """..."""
        self._axis_update("_redistribution_matrix_input_axis_edges", 
                          input_axis, 
                          "Redistribution matrix input axis")

    def _set_rmf_output_axis(self, output_axis):
        """..."""
        self._axis_update("_redistribution_matrix_output_axis_edges", 
                          output_axis, 
                          "Redistribution matrix output axis")
            
    def get_srm_info(self):
        """..."""
        if self._has_srm:
            return {"srm_input_axis":self._spectral_response_matrix_input_axis_edges,
                    "srm_output_axis":self._spectral_response_matrix_output_axis_edges,
                    "srm":self._spectral_response_matrix}
        
        warnings.warn("Missing SRM information.")

    def set_srm_info(self, srm=None, output_axis_edges=None, input_axis_edges=None):
        """..."""
        warnings.warn(
        """
        Be careful updating the SRM directly. If you really have to do 
        something it is recommended to edit the ARF or RMF separately 
        then construct a new SRM otherwise the maths will likely not
        work in your favour.
        """
        )

        if srm is not None:
            self._set_srm(srm)
        if output_axis_edges is not None:
            warnings.warn("""
            Editing the output axis will affect the SRM's compatibility 
            with the PHA, spectrum data.
            """
            )
            self._set_srm_output_axis(output_axis_edges)
        if input_axis_edges is not None:
            warnings.warn("""
            OK, you really should not be editing the input axis here. Go 
            update the ARF and RMF separately then construct a new SRM 
            but, hey, what do I know!
            """
            )
            self._set_srm_input_axis(input_axis_edges)

        self._post_srm_update_checks()

        if self._has_srm:
            warnings.warn("Overwritten SRM information.")
        self._has_srm = True

    def _set_srm(self, srm):
        """..."""
        self._spectral_response_matrix = srm

        self.standard_unit_check("Spectral response matrix", 
                                 self._spectral_response_matrix.unit, 
                                 "srm")

    def _axis_update(self, att, new_val, desc):
        """..."""
        if not self._has_unit(att, new_val, desc, "energy"):
            return
        self.__dict__[att] = new_val

    def _has_unit(self, att, new_val, desc, standard_unit):
        if not isinstance(new_val, u.Quantity):
            warnings.warn(f"""
            Nom, nom, nom, give me yummy units on `{att}` update 
            otherwise I'm not changing a thing.
            """
            )
            return False

        self.standard_unit_check(desc, 
                                 self.__dict__[att].unit, 
                                 standard_unit)
        return True

    def _set_srm_input_axis(self, input_axis):
        """..."""
        self._axis_update("_spectral_response_matrix_input_axis_edges", 
                          input_axis, 
                          "Spectral response matrix input axis")

    def _set_srm_output_axis(self, output_axis):
        """..."""
        self._axis_update("_spectral_response_matrix_output_axis_edges", 
                          output_axis, 
                          "Spectral response matrix output axis")

    def _post_srm_update_checks(self):
        # check the shape of the SRM at least
        srm_shape = np.shape(self._spectral_response_matrix)
        if srm_shape[0]!=np.shape(self._spectral_response_matrix_input_axis_edges)[0]:
            warnings.warn("Input photon axis length of SRM does not match input axis values.")
        if srm_shape[1]!=np.shape(self._spectral_response_matrix_output_axis_edges)[0]:
            warnings.warn("Output count axis length of SRM does not match output axis values.")


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

        (chan, e_min, e_max), (e_lo_rmf, e_hi_rmf, ngrp, fchan, nchan, matrix) = get_response_info(*read_heasarc_rmf(rmf_file))
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
