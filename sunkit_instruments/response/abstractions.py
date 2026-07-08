"""This module defines abstractions for computing instrument response."""
import abc
import copy
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

    A channel instance is a snapshot of the instrument at one time:
    time-dependent degradation enters exclusively through `at`, and every
    consumer-facing method (`effective_area`, `wavelength_response`) takes
    no time argument. Subclasses implement the time dependence in
    `degradation`.

    For all methods and properties defined here, see the
    topic guide on instrument response for more information.
    """

    _obstime = None

    def at(self, obstime):
        """
        A copy of this channel evaluated at ``obstime``.

        The single way time enters a response calculation: bind the time
        first, then pass the returned channel wherever a channel is
        accepted::

            get_temperature_response(channel.at("2020-01-01"), model)

        Parameters
        ----------
        obstime: any format parsed by `sunpy.time.parse_time`
        """
        bound = copy.copy(self)
        bound._obstime = obstime
        return bound

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
        time bound with `at`. An unbound channel (including ``at(None)``)
        never evaluates degradation — it is the pristine instrument.
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
