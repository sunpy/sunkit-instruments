import warnings

import numpy as np
from ndcube import NDMeta
from sunkit_spex.spectrum.spectrum import SpectralAxis, Spectrum
from sunkit_spex.spectrum.uncertainty import PoissonUncertainty

import astropy.units as u

from sunkit_instruments.nustar.io import (
    read_heasarc_arf,
    read_heasarc_rmf,
    read_nustar_pha,
)
from sunkit_instruments.nustar.spectrum import (
    col2arr,
    get_effective_area_info,
    get_observable_info,
    get_response_info,
    make_srm,
    vrmf2arr,
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
            Default: None
    """

    def __init__(self, pha_file:str|None=None, arf_file:str|None=None, rmf_file:str|None=None, **kwargs):
        """Define standard units, a construction string for repr, and construct."""
        self._construction_string = (
            f"{str(self.__class__)}(pha_file={pha_file},arf_file={arf_file},rmf_file={rmf_file})"
        )
        self.define_standard_units()
        # RMF contains channel->count energy information so do this first if we can
        self.construct_rmf(f_rmf=rmf_file)
        self.construct_pha(pha_file)
        self.construct_arf(f_arf=arf_file)
        # the SRM construction needs the RMF and ARF to be done
        self.construct_srm()
        self._spec_obj_inputs = {"distance":1.0<<u.AU} | kwargs
        self.construct_spec_obj()

    def define_standard_units(self):
        """Method to define standard units for everything stored in the class."""
        self._standard_units = {"channel_number":(u.dimensionless_unscaled),
                                "energy":(u.keV),
                                "ct_spec":(u.ct),
                                "eff_exp/lvt":(u.s),
                                "eff_area":(u.cm**2),
                                "rdm":(u.ct * u.ph**-1),
                                "srm":(u.ct * u.ph**-1 * u.cm**2)}

    def construct_spec_obj(self, **kwargs):
        """Will construct a spectrum object for Sunkit-spex.

        Creates the ``spectrum_object`` attribute. The `kwargs` are
        passed to the returned object's meta.
        """
        self._has_spec_obj = False
        if not self._has_pha:
            warnings.warn("No PHA information available to create spectrum object.")
            return

        pha_info = self.get_pha_info()
        counts = pha_info["spectrum_counts"]
        counts_uncertainity = np.sqrt(counts.value)<<counts.unit
        count_axis = pha_info["spectrum_counts_axis"]
        exp_time = pha_info["effective_exposure"]

        counts_uncertainity_pu = PoissonUncertainty(counts_uncertainity)
        count_axis_flat = np.append(count_axis[:,0], count_axis[-1,1])
        counts_spectral_axis = SpectralAxis(count_axis_flat, bin_specification="edges")

        meta = NDMeta()
        meta.add("exposure_time", exp_time)

        if self._has_srm:
            srm_info = self.get_srm_info()
            meta.add("srm", srm_info["srm"])
            meta.add("ph_axis", srm_info["srm_input_axis"])

        self._spec_obj_inputs |= kwargs

        for (k, v) in self._spec_obj_inputs.items():
            meta.add(k, v)

        self.spectrum_object = Spectrum(
            data=counts, uncertainty=counts_uncertainity_pu, spectral_axis=counts_spectral_axis, meta=meta
        )
        self._has_spec_obj = True

    def get_spec_obj(self):
        """Return the spectrum object if it exists, else create one."""
        if self._has_spec_obj:
            return self.spectrum_object

        self.construct_spec_obj()
        if hasattr(self, "spectrum_object"):
            return self.spectrum_object

    def _get_observable_info(self, f_pha:str):
        """Separate function for easy testing."""
        return get_observable_info(*read_nustar_pha(f_pha))

    def construct_pha(self, f_pha:str|None=None):
        """Read a `.pha` file and store all the spectral information."""
        if f_pha is None:
            warnings.warn(f"File `{f_pha}` has not been given.")
            self._spectrum_channel_number = None
            self._spectrum_counts = None
            self._effective_exposure = None
            self._has_pha = False
            return

        self._spectrum_channel_number, self._spectrum_counts, self._effective_exposure = self._get_observable_info(f_pha)
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
            warnings.warn("No RMF information exists so defaulting to standard NuSTAR energy binning for the spectrum.")
            self.set_standard_spectrum_axis()

    def set_standard_spectrum_axis(self):
        """Defines the usual, native NuSTAR energy binning."""
        _standard_energy_start = 1.6 << u.keV
        _standard_num_of_channels = 4096
        _standard_energy_binning = 0.04 << u.keV
        _standard_max_energy = _standard_num_of_channels*_standard_energy_binning + _standard_energy_start
        e_lo = np.arange(_standard_energy_start.value, _standard_max_energy.value, _standard_energy_binning.value) << _standard_energy_start.unit
        e_hi = e_lo + _standard_energy_binning
        self._spectrum_axis_edges = np.hstack((e_lo[:,None], e_hi[:,None]))

    def _get_effective_area_info(self, f_arf:str):
        """Separate function for easy testing."""
        return get_effective_area_info(read_heasarc_arf(f_arf))

    def construct_arf(self, f_arf:str|None=None):
        """Read an `.arf` file and store all the response information."""
        if f_arf is None:
            warnings.warn(f"File `{f_arf}` has not been given.")
            self._effective_area_axis_edges = None
            self._effective_area = None
            self._has_arf = False
            return

        e_lo_arf, e_hi_arf, self._effective_area = self._get_effective_area_info(f_arf)
        self._effective_area_axis_edges = np.hstack((e_lo_arf[:,None], e_hi_arf[:,None]))
        self._has_arf = True

        self.standard_unit_check("Effective area array",
                                 self._effective_area.unit,
                                 "eff_area")
        self.standard_unit_check("Effective area axis",
                                 self._effective_area_axis_edges.unit,
                                 "energy")

    def construct_rmf(self, f_rmf:str|None=None):
        """Read a `.rmf` file and store all the matrix information."""
        if f_rmf is None:
            warnings.warn(f"File `{f_rmf}` has not been given.")
            self._redistribution_matrix_aux_info = None
            self._redistribution_matrix_input_axis_edges = None
            self._redistribution_matrix_ouput_channel_number = None
            self._redistribution_matrix_output_axis_edges = None
            self._redistribution_matrix = None
            self._has_rmf = False
            return

        (chan, e_min, e_max), (e_lo_rmf, e_hi_rmf, ngrp, fchan, nchan, matrix, redist_m) = self._load_rmf(f_rmf)

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
        """Create a spectral response matrix from the ARF and RMF information."""
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
        """Return the PHA, spectrum information, if available."""
        if self._has_pha:
            return {"spectrum_counts_axis":self._spectrum_axis_edges,
                    "spectrum_counts":self._spectrum_counts,
                    "effective_exposure":self._effective_exposure}

        warnings.warn("Missing PHA information.")

    def set_pha_info(self, spectrum_counts:u.Quantity|None=None, spectrum_counts_axis:u.Quantity|None=None, effective_exposure:u.Quantity|None=None, pha_file:str|None=None):
        """Allows parts of the PHA, spectrum information to be updated.

        Parameters
        ----------
        spectrum_counts : `~astropy.units.Quantity`
            The observable of the NuSTAR spectrum. Normally, this is in
            counts.
            Default: None

        spectrum_counts_axis : `~astropy.units.Quantity`
            The bin edges on which `spectrum_counts` is measured. Normally,
            this is in keV.
            Default: None

        effective_exposure : `~astropy.units.Quantity`
            The effective exposure for the spectrum. Normally, this is in
            seconds.
            Default: None

        pha_file : `str`
            A new `.pha` file to update the spectrum information. This
            takes priority over all other arguments.
            Default: None
        """
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
            if effective_exposure is not None:
                self._set_pha_effective_exposure(effective_exposure)
        else:
            self.construct_pha()

        self._post_pha_update_checks()

        if self._has_pha:
            warnings.warn("Overwritten PHA information.")
        self._has_pha = True

    def _post_pha_update_checks(self):
        """A function to check updated PHA information."""
        # check the shapes match at least
        pha_shape = np.shape(self._spectrum_counts)
        if pha_shape[0]!=np.shape(self.spectrum_axis_edges)[0]:
            warnings.warn("Count axis length of PHA does not match PHA length.")

    def _set_pha(self, pha:u.Quantity):
        """Sets a new spectrum array and checks units."""
        self._spectrum_counts = pha

        self.standard_unit_check("Count spectrum array",
                                 self._spectrum_counts.unit,
                                 "ct_spec")

    def _set_pha_axis(self, axis:u.Quantity):
        """Updates the spectrum bin edges."""
        self._axis_update("_spectrum_axis_edges",
                          axis,
                          "Spectrum energy axis")

    def _set_pha_effective_exposure(self, time:u.Quantity):
        """Sets a new effective exposure value and checks units."""
        self._effective_exposure = time

        self.standard_unit_check("Effective exposure/livetime value",
                                 self._effective_exposure.unit,
                                 "eff_exp/lvt")

    def get_arf_info(self):
        """Return the ARF, effective area information, if available."""
        if self._has_arf:
            return {"effective_area_axis":self._effective_area_axis_edges,
                    "effective_area":self._effective_area}

        warnings.warn("Missing ARF information.")

    def set_arf_info(self, effective_area:u.Quantity|None=None, effective_area_axis:u.Quantity|None=None, arf_file:str|None=None):
        """Allows parts of the ARF, effective area information to be updated.

        Parameters
        ----------
        effective_area : `~astropy.units.Quantity`
            The effective area array. Normally, this is in cm.
            Default: None

        effective_area_axis : `~astropy.units.Quantity`
            The bin edges on which `effective_area` is measured. Normally,
            this is in keV.
            Default: None

        arf_file : `str`
            A new `.arf` file to update the effective area information.
            This takes priority over all other arguments.
            Default: None
        """
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
        """A function to check updated ARF information."""
        # check the shapes match at least
        arf_shape = np.shape(self._effective_area)
        if arf_shape[0]!=np.shape(self._effective_area_axis_edges)[0]:
            warnings.warn("Photon axis length of ARF does not match ARF length.")

    def _set_arf(self, arf:u.Quantity):
        """Sets a new effective area array and checks units."""
        self._effective_area = arf

        self.standard_unit_check("Effective area array",
                                 self._effective_area.unit,
                                 "eff_area")

    def _set_arf_axis(self, axis:u.Quantity):
        """Updates the effective area bin edges."""
        self._axis_update("_effective_area_axis_edges",
                          axis,
                          "Effective area axis")

    def get_rmf_info(self, include_auxilliray_info:bool=False):
        """Return the RMF, redistribution matrix information, if available.

        Parameters
        ----------
        include_auxilliray_info : `bool`
            Updates the returning dictionary with the lower level data
            used to create the redistribution matrix.
            Default: False
        """
        if self._has_rmf:
            info = {"rmf_input_axis":self._redistribution_matrix_input_axis_edges,
                    "rmf_output_axis":self._redistribution_matrix_output_axis_edges,
                    "rmf":self._redistribution_matrix}

            if include_auxilliray_info:
                return info | self._redistribution_matrix_aux_info

            return info

        warnings.warn("Missing RMF information.")

    def set_rmf_info(self, rmf:u.Quantity|None=None, output_axis_edges:u.Quantity|None=None, input_axis_edges:u.Quantity|None=None, rmf_file:str|None=None):
        """Allows parts of the RMF, redistribution matrix information to be updated.

        Parameters
        ----------
        rmf : `~astropy.units.Quantity`
            The redistribution matrix. Normally, this is in counts/photon.
            Default: None

        output_axis_edges : `~astropy.units.Quantity`
            The bin edges on which `rmf` columns are defined. Normally,
            this is in keV.
            Default: None

        input_axis_edges : `~astropy.units.Quantity`
            The bin edges on which `rmf` rows are defined. Normally,
            this is in keV.
            Default: None

        rmf_file : `str`
            A new `.rmf` file to update the effective area information.
            This takes priority over all other arguments.
            Default: None
        """
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
        """A function to check updated RMF information."""
        # check the shapes match at least
        rmf_shape = np.shape(self._redistribution_matrix)
        if rmf_shape[0]!=np.shape(self._redistribution_matrix_input_axis_edges)[0]:
            warnings.warn("Input photon axis length of RMF does not match input axis values.")
        if rmf_shape[1]!=np.shape(self._redistribution_matrix_output_axis_edges)[0]:
            warnings.warn("Output count axis length of RMF does not match output axis values.")

    def _set_rmf(self, rmf:u.Quantity):
        """Sets a new redistribution matrix and checks units."""
        self._redistribution_matrix = rmf

        self.standard_unit_check("Redistribution matrix",
                                 self._redistribution_matrix.unit,
                                 "rmf")

    def _set_rmf_input_axis(self, input_axis:u.Quantity):
        """Updates the redistribution matrix input (rows) bin edges."""
        self._axis_update("_redistribution_matrix_input_axis_edges",
                          input_axis,
                          "Redistribution matrix input axis")

    def _set_rmf_output_axis(self, output_axis:u.Quantity):
        """Updates the redistribution matrix output (columns) bin edges."""
        self._axis_update("_redistribution_matrix_output_axis_edges",
                          output_axis,
                          "Redistribution matrix output axis")

    def get_srm_info(self):
        """Return the spectral response information, if available."""
        if self._has_srm:
            return {"srm_input_axis":self._spectral_response_matrix_input_axis_edges,
                    "srm_output_axis":self._spectral_response_matrix_output_axis_edges,
                    "srm":self._spectral_response_matrix}

        warnings.warn("Missing SRM information.")

    def set_srm_info(self, srm:u.Quantity|None=None, output_axis_edges:u.Quantity|None=None, input_axis_edges:u.Quantity|None=None):
        """Allows parts of the SRM, spectral response matrix information to be updated.

        Parameters
        ----------
        srm : `~astropy.units.Quantity`
            The spectral response matrix. Normally, this is in
            (counts * cm^2)/photon.
            Default: None

        output_axis_edges : `~astropy.units.Quantity`
            The bin edges on which `arm` columns are defined. Normally,
            this is in keV.
            Default: None

        input_axis_edges : `~astropy.units.Quantity`
            The bin edges on which `arm` rows are defined. Normally,
            this is in keV.
            Default: None
        """
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

    def _set_srm(self, srm:u.Quantity):
        """Sets a new spectrum response matrix and checks units."""
        self._spectral_response_matrix = srm

        self.standard_unit_check("Spectral response matrix",
                                 self._spectral_response_matrix.unit,
                                 "srm")

    def _axis_update(self, att:str, new_val:u.Quantity, desc:str):
        """General method to update axes with units of energy."""
        if not self._has_unit(att, new_val, desc, "energy"):
            return
        self.__dict__[att] = new_val

    def _has_unit(self, att:str, new_val:u.Quantity, desc:str, standard_unit:str):
        """General method to check if a new value given has a unit."""
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
        """Updates the spectral response matrix input (rows) bin edges."""
        self._axis_update("_spectral_response_matrix_input_axis_edges",
                          input_axis,
                          "Spectral response matrix input axis")

    def _set_srm_output_axis(self, output_axis):
        """Updates the spectral response matrix output (columns) bin edges."""
        self._axis_update("_spectral_response_matrix_output_axis_edges",
                          output_axis,
                          "Spectral response matrix output axis")

    def _post_srm_update_checks(self):
        """A function to check updated SRM information."""
        # check the shape of the SRM at least
        srm_shape = np.shape(self._spectral_response_matrix)
        if srm_shape[0]!=np.shape(self._spectral_response_matrix_input_axis_edges)[0]:
            warnings.warn("Input photon axis length of SRM does not match input axis values.")
        if srm_shape[1]!=np.shape(self._spectral_response_matrix_output_axis_edges)[0]:
            warnings.warn("Output count axis length of SRM does not match output axis values.")

    def get_count_energy_bins(self):
        """Return the count bin edges (SRM and RMF's output), if available."""
        if not self._has_rmf:
            warnings.warn("Missing count energy bins as no RMF information is available. Returning default spectrum axes instead.")
            return self._spectrum_axis_edges

        if not np.allclose(self._spectrum_axis_edges, self._redistribution_matrix_output_axis_edges):
            warnings.warn("Spectrum and RMF count bin edge information is different. Returning dictionary of both.")
            return {"spec-count-bin-edges":self._spectrum_axis_edges,
                    "rmf-count-bin-edges":self._redistribution_matrix_input_axis_edges}

        return self._redistribution_matrix_output_axis_edges

    def get_photon_energy_bins(self):
        """Return the photon bin edges (SRM and RMF's input), if available."""
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

    def standard_unit_check(self, name:str, unit:u.core.PrefixUnit|u.core.CompositeUnit, key:str):
        """Checks units against the ``_standard_units`` attribute."""
        if unit!=self._standard_units[key]:
            warnings.warn(f"{name} units are not standard {self._standard_units[key]}, but in {unit}.")

    def _get_response_info(self, rmf_file:str):
        """Separate function for easy testing."""
        return get_response_info(*read_heasarc_rmf(rmf_file))

    def _load_rmf(self, rmf_file:str):
        """Extracts all information, mainly the redistribution matrix
        ([counts/photon]) from a given RMF file.

        Parameters
        ----------
        rmf_file : `str`
            The file path and name of the RMF file.

        Returns
        -------
        Two tuples:

        (1) The channel number and their associate lower and
        upper energy bounds.

        (2) The lower/higher photon bin edges (e_lo_rmf, e_hi_rmf), the
        number of counts channels activated by each photon channel (ngrp),
        starting indices of the count channel groups (fchan), number counts
        channels from each starting index (nchan), the corresponding
        counts/photon value for each count and photon entry (matrix), and
        the redistribution matrix (redist_m: with rows of photon channels,
        columns of counts channels, and in the units of counts/photon).
        """

        (chan, e_min, e_max), (e_lo_rmf, e_hi_rmf, ngrp, fchan, nchan, matrix) = self._get_response_info(rmf_file)
        fchan_array = col2arr(fchan)
        nchan_array = col2arr(nchan)
        redist_m = vrmf2arr(
            data=matrix, n_grp_list=ngrp, f_chan_array=fchan_array, n_chan_array=nchan_array
        )

        return (chan, e_min, e_max), (e_lo_rmf, e_hi_rmf, ngrp, fchan, nchan, matrix, redist_m)

    def rebin_rmf_info(self, new_input_axis_edges:u.Quantity|None=None, new_output_axis_edges:u.Quantity|None=None):
        """Function to rebin the axes of the RMF, redistribution matrix.

        The PHA, ARF, RMF, and SRM should all be updated together so
        please used: ``rebin_info`` instead.

        Parameters
        ----------
        new_input_axis_edges : `~astropy.units.Quantity`
            The new bin edges for the input (rows) RMF axis. This axis is
            shared with the ARF axis.
            Default: None

        new_output_axis_edges : `~astropy.units.Quantity`
            The new bin edges for the output (columns) RMF axis. This axis
            is shared with the PHA axis.
            Default: None

        Returns
        -------
        : `tuple[~astropy.units.Quantity]`
            New input axis edges, output axis edges, and the new RMF.
        """
        warnings.warn("Only rebinning the RMF information is ill-advised.")

        rmf_info = self.get_rmf_info()

        # checked with ftrbnrmf
        new_rmf = rebin_rmf(
            matrix=rmf_info["rmf"],
            old_output_bins=rmf_info["rmf_output_axis"],
            new_output_bins=new_output_axis_edges,
            old_input_bins=rmf_info["rmf_input_axis"],
            new_input_bins=new_input_axis_edges
        )
        new_input_axis_edges = new_input_axis_edges if new_input_axis_edges is not None else rmf_info["rmf_input_axis"]
        new_output_axis_edges = new_output_axis_edges if new_output_axis_edges is not None else rmf_info["rmf_output_axis"]
        return new_input_axis_edges, new_output_axis_edges, new_rmf

    def rebin_arf_info(self, new_axis_edges:u.Quantity):
        """Function to rebin the axes of the ARF, effective area information.

        The PHA, ARF, RMF, and SRM should all be updated together so
        please used: ``rebin_info`` instead.

        Parameters
        ----------
        new_axis_edges : `~astropy.units.Quantity`
            The new bin edges for the ARF axis. This axis is shared with
            the RMF/SRM input axis.

        Returns
        -------
        : `tuple[~astropy.units.Quantity]`
            New axis edges and the new ARF.
        """
        warnings.warn("Only rebinning the ARF information is ill-advised.")

        arf_info = self.get_arf_info()

        new_arf = regroup_any_array(
            data=arf_info["effective_area"],
            old_bins=arf_info["effective_area_axis"],
            new_bins=new_axis_edges,
            combine_by="mean"
        )
        return new_axis_edges, new_arf

    def rebin_pha_info(self, new_axis_edges:u.Quantity):
        """Function to rebin the axes of the PHA, spectrum information.

        The PHA, ARF, RMF, and SRM should all be updated together so
        please used: ``rebin_info`` instead.

        Parameters
        ----------
        new_axis_edges : `~astropy.units.Quantity`
            The new bin edges for the PHA axis. This axis is shared with
            the RMF/SRM output axis.

        Returns
        -------
        : `tuple[~astropy.units.Quantity]`
            New axis edges and the new PHA.
        """
        warnings.warn("Only rebinning the PHA information is ill-advised.")

        pha_info = self.get_pha_info()

        new_pha = regroup_any_array(
            data=pha_info["spectrum_counts"],
            old_bins=pha_info["spectrum_counts_axis"],
            new_bins=new_axis_edges,
            combine_by="sum"
        )
        return new_axis_edges, new_pha

    def rebin_srm_info(self, new_input_axis_edges:u.Quantity|None=None, new_output_axis_edges:u.Quantity|None=None):
        """Function to rebin the axes of the SRM, spectral response matrix.

        The PHA, ARF, RMF, and SRM should all be updated together so
        please used: ``rebin_info`` instead.

        Parameters
        ----------
        new_input_axis_edges : `~astropy.units.Quantity`
            The new bin edges for the input (rows) SRM axis. This axis is
            shared with the ARF axis and input RMF axis.
            Default: None

        new_output_axis_edges : `~astropy.units.Quantity`
            The new bin edges for the output (columns) SRM axis. This axis
            is shared with the PHA axis and output RMF axis.
            Default: None

        Returns
        -------
        : `tuple[~astropy.units.Quantity]`
            New input axis edges, output axis edges, the new ARF, the new
            RMF, and the new SRM.
        """
        warnings.warn("Only rebinning the SRM information, and not PHA as well, is ill-advised.")
        with warnings.catch_warnings(action="ignore"):
            new_input, new_output, new_rmf = self.rebin_rmf_info(new_input_axis_edges=new_input_axis_edges,
                                                new_output_axis_edges=new_output_axis_edges)
            _, new_arf = self.rebin_arf_info(new_axis_edges=new_input_axis_edges)

        new_srm = make_srm(rmf_matrix=new_rmf, arf_array=new_arf)

        new_input_axis_edges = new_input_axis_edges if new_input_axis_edges is not None else new_input
        new_output_axis_edges = new_output_axis_edges if new_output_axis_edges is not None else new_output

        return new_input_axis_edges, new_output_axis_edges, new_arf, new_rmf, new_srm

    def rebin_info(self, new_input_axis_edges:u.Quantity|None=None, new_output_axis_edges:u.Quantity|None=None):
        """Function to rebin the axes of the ARF, RMF (so SRM as well),
        and PHA information.

        Parameters
        ----------
        new_input_axis_edges : `~astropy.units.Quantity`
            The new bin edges for the input (rows) SRM axis. This axis is
            shared with the ARF axis and input RMF axis.
            Default: None

        new_output_axis_edges : `~astropy.units.Quantity`
            The new bin edges for the output (columns) SRM axis. This axis
            is shared with the PHA axis and output RMF axis.
            Default: None

        Returns
        -------
        : `tuple[~astropy.units.Quantity]`
            ...
        """
        with warnings.catch_warnings(action="ignore"):
            if new_output_axis_edges is not None:
                new_axis_edges, new_pha = self.rebin_pha_info(new_axis_edges=new_output_axis_edges)
            new_input_axis_edges, new_output_axis_edges, new_arf, new_rmf, new_srm = self.rebin_srm_info(
                new_input_axis_edges=new_input_axis_edges,
                new_output_axis_edges=new_output_axis_edges
                )

            self.set_pha_info(spectrum_counts=new_pha,
                              spectrum_counts_axis=new_axis_edges)
            self.set_arf_info(effective_area=new_arf,
                              effective_area_axis=new_input_axis_edges)
            self.set_rmf_info(rmf=new_rmf,
                              output_axis_edges=new_output_axis_edges,
                              input_axis_edges=new_input_axis_edges)
            self.set_srm_info(srm=new_srm,
                              output_axis_edges=new_output_axis_edges,
                              input_axis_edges=new_input_axis_edges)

    def __repr__(self):
        """String representation of the class."""
        return self._construction_string
