"""This module defines abstractions for computing instrument response."""
import abc

import numpy as np
from scipy.interpolate import interp1d

import astropy.units as u

__all__ = ["AbstractChannel"]


class AbstractChannel(abc.ABC):
    """
    An abstract base class for defining instrument channels.

    For all methods and properties defined here, see the
    `topic guide on instrument response <>`_ for more information.
    """

    @u.quantity_input
    def temperature_response(
        self, source_spectra, obstime=None
    ) -> u.cm**5 * u.DN / (u.pix * u.s):
        """
        Temperature response function for a given source spectra.

        The temperature response function describes the sensitivity of an imaging
        instrument as a function of temperature. The temperature response is
        calculated by integrating the source spectra over the wavelength dimension,
        weighted by the wavelength response of the instrument.

        Parameters
        ----------
        source_spectra: `SourceSpectra`
        obstime: any format parsed by `~sunpy.time.parse_time`, optional
        """
        # TODO: refactor all of this to take advantage of xarray interpolation
        # and summation
        wave_response = self.wavelength_response(obstime=obstime)
        f_response = interp1d(
            self.wavelength.to_value("AA"),
            wave_response.to_value(),
            axis=0,
            bounds_error=False,
            fill_value=0.0,
        )
        wave_response_interp = u.Quantity(
            f_response(source_spectra.wavelength.to_value("AA")), wave_response.unit
        )
        temperature_response = (
            source_spectra.data
            * wave_response_interp
            * np.gradient(source_spectra.wavelength)
        ).sum(axis=1)
        return temperature_response

    @u.quantity_input
    def wavelength_response(
        self, obstime=None
    ) -> u.cm**2 * u.DN * u.steradian / (u.photon * u.pixel):
        """
        Instrument response as a function of wavelength

        The wavelength response is the effective area with
        the conversion factors from photons to DN and steradians
        to pixels.

        Parameters
        ----------
        obstime: any format parsed by `~sunpy.time.parse_time`, optional
            If specified, this is used to compute the time-dependent
            instrument degradation.
        """
        area_eff = self.effective_area(obstime=obstime)
        return (
            area_eff
            * self.energy_per_photon
            * self.pixel_solid_angle
            * self.camera_gain
            / self.energy_per_electron
        )

    @u.quantity_input
    def effective_area(self, obstime=None) -> u.cm**2:
        """
        Effective area as a function of wavelength.

        The effective area is the geometrical collecting area
        weighted by the mirror reflectance, filter transmittance,
        quantum efficiency, and instrument degradation.

        Parameters
        ----------
        obstime: any format parsed by `~sunpy.time.parse_time`, optional
            If specified, this is used to compute the time-dependent
            instrument degradation.
        """
        return (
            self.geometrical_area
            * self.mirror_reflectance
            * self.filter_transmittance
            * self.effective_quantum_efficiency
            * self.degradation(obstime=obstime)
        )

    @property
    @u.quantity_input
    def energy_per_photon(self) -> u.eV / u.photon:
        return self.wavelength.to("eV", equivalencies=u.spectral()) / u.photon

    @abc.abstractmethod
    @u.quantity_input
    def degradation(self, obstime=None): ...

    @property
    @abc.abstractmethod
    @u.quantity_input
    def geometrical_area(self) -> u.cm**2: ...

    @property
    @abc.abstractmethod
    @u.quantity_input
    def mirror_reflectance(self) -> u.dimensionless_unscaled: ...

    @property
    @abc.abstractmethod
    @u.quantity_input
    def filter_transmittance(self) -> u.dimensionless_unscaled: ...

    @property
    @abc.abstractmethod
    @u.quantity_input
    def effective_quantum_efficiency(self) -> u.dimensionless_unscaled: ...

    @property
    @abc.abstractmethod
    @u.quantity_input
    def camera_gain(self) -> u.DN / u.electron: ...

    @property
    @abc.abstractmethod
    @u.quantity_input
    def energy_per_electron(self) -> u.eV / u.electron: ...

    @property
    @abc.abstractmethod
    @u.quantity_input
    def pixel_solid_angle(self) -> u.steradian / u.pixel: ...

    @property
    @abc.abstractmethod
    @u.quantity_input
    def wavelength(self) -> u.Angstrom: ...
