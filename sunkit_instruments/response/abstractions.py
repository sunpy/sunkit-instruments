"""This module defines abstractions for computing instrument response."""
import abc
from typing import Protocol, runtime_checkable

import astropy.units as u

__all__ = ["AbstractChannel", "EmissionModel", "LineEmissionModel"]


@runtime_checkable
class EmissionModel(Protocol):
    """
    Protocol for emission models usable with
    `~sunkit_instruments.response.get_temperature_response`.

    An emission model bundles the atomic physics: it knows how a plasma
    emits as a function of temperature and can fold that emission through
    an instrument channel. Any object providing this interface conforms —
    no inheritance from this class is required. In particular
    `synthesizAR.atomic.EmissionModel` satisfies it.

    Notes
    -----
    Implementations should document whether the emissivity used by
    `calculate_temperature_response` already includes the (equilibrium)
    ionization fraction, since some models (e.g. synthesizAR) defer it
    to support time-dependent ionization.
    """

    @property
    def temperature(self) -> u.Quantity:
        """Temperature grid of the model."""

    def calculate_temperature_response(self, channel) -> u.Quantity:
        """
        Temperature response of ``channel`` for this model's emission.

        ``channel`` is compatible with
        `~sunkit_instruments.response.abstractions.AbstractChannel`; its
        ``wavelength`` and ``wavelength_response()`` are the interface
        implementations should rely on.
        """


@runtime_checkable
class LineEmissionModel(Protocol):
    """
    Protocol for emission models that expose per-transition line emissivities.

    This is the surface consumed by
    `~sunkit_instruments.response.line_list_from_emission_model` to build a
    line list for spectrally-resolved response functions
    (`~sunkit_instruments.response.create_response_function`). The model is
    an iterable of ion objects (each providing ``atomic_number``,
    ``ion_name_roman`` and ``ionization_fraction``), matching
    `synthesizAR.atomic.EmissionModel` and `fiasco.IonCollection`.
    """

    @property
    def temperature(self) -> u.Quantity:
        """Temperature grid of the model."""

    @property
    def density(self) -> u.Quantity:
        """Density grid of the model."""

    def __iter__(self):
        """Iterate over the ions in the model."""

    def get_line_emissivity(self, ion, transition=None):
        """
        Per-transition line emissivity for one ion.

        Returns a ``(wavelength, emissivity)`` pair where ``emissivity`` has
        shape ``(temperature, density, transition)`` in photon units
        (``cm3 ph s-1``) and excludes the ionization fraction, matching the
        synthesizAR convention.
        """


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
