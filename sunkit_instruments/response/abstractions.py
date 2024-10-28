"""This module defines abstractions for computing instrument response."""
import abc

import astropy.units as u

__all__ = ["AbstractChannel"]


class AbstractChannel(abc.ABC):
    """
    An abstract base class for defining instrument channels.

    For all methods and properties defined here, see the
    topic guide on instrument response for more information.
    """

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
        obstime: any format parsed by `sunpy.time.parse_time`, optional
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
    def degradation(self, obstime=None) -> u.dimensionless_unscaled: ...

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
