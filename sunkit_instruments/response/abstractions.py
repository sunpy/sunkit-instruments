"""This module defines abstractions for computing instrument response."""
import abc
import copy

import astropy.units as u

from sunpy.time import parse_time

__all__ = ["AbstractChannel"]


class AbstractChannel(abc.ABC):
    """
    An abstract base class for defining instrument channels.

    A channel instance can be a snapshot of the instrument at one time:
    time-dependent degradation enters exclusively through the ``at`` method.
    Subclasses implement the time dependence in ``degradation``.

    For all methods and properties defined here, see the
    topic guide on instrument response for more information.
    """

    _obstime = None

    def at(self, obstime):
        """
        A copy of this channel evaluated at ``obstime``.

        The single way time enters any calculation.
        Bind the time first, then pass the returned channel wherever a channel is
        accepted:

        .. code-block:: python

            get_temperature_response(channel.at("2020-01-01"), spectra)

        Parameters
        ----------
        obstime: any format parsed by `sunpy.time.parse_time`
        """
        bound = copy.copy(self)
        bound._obstime = parse_time(obstime) if obstime is not None else None
        return bound

    def __repr__(self):
        obstime = "Non-degraded" if self._obstime is None else f"Degraded at {self._obstime.isot}"
        return f"{type(self).__name__}({obstime})"

    @u.quantity_input
    def wavelength_response(self) -> u.cm**2 * u.DN * u.steradian / (u.photon * u.pixel):
        """
        Instrument response as a function of wavelength

        The wavelength response is the effective area with
        the conversion factors from photons to DN and steradians
        to pixels. For time-dependent degradation, bind the time first
        with `at`.
        """
        return (
            self.effective_area()
            * self.energy_per_photon
            * self.pixel_solid_angle
            * self.camera_gain
            / self.energy_per_electron
        )

    @u.quantity_input
    def effective_area(self) -> u.cm**2:
        """
        Effective area as a function of wavelength.

        The effective area is the geometrical collecting area
        weighted by the mirror reflectance, filter transmittance,
        quantum efficiency, and instrument degradation evaluated at the
        time bound with ``at``. An unbound channel never evaluates degradation.
        """
        area = (
            self.geometrical_area
            * self.mirror_reflectance
            * self.filter_transmittance
            * self.effective_quantum_efficiency
        )
        if self._obstime is None:
            return area
        return area * self.degradation(obstime=self._obstime)

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
