"""This module defines abstractions related to instrument channels."""
import abc

import numpy as np
import xarray
from scipy.interpolate import interp1d

import astropy.units as u

__all__ = ["AbstractChannel", "SourceSpectra"]


class SourceSpectra:
    """
    Source spectra as a function of temperature and wavelength.

    The source spectra describes how a plasma emits (under the optically-thin
    assumption) as a function of both temperature and wavelength. This source
    spectra is typically computed using a database like CHIANTI by summing the
    emission spectra of many ions as well as the continuum emission. For more
    information, see the `topic guide on instrument response <>`_.

    Parameters
    ----------
    temperature: u.Quantity
        1D array describing the variation along the temperature axis
    wavelength: u.Quantity
        1D array describing the variation along the wavelength axis
    spectra: u.Quantity
        Source spectra as a 2D array. The first axis should correspond to temperature and the second axis should correspond to wavelength.
    density: u.Quantity, optional
        1D array describing the variation in density along the temperature axis. It is assumed
        that temperature and density are dependent.
    """

    @u.quantity_input
    def __init__(
        self,
        temperature: u.K,
        wavelength: u.Angstrom,
        spectra: u.photon * u.cm**3 / (u.s * u.AA * u.sr),
        density: u.cm ** (-3) = None,
        **kwargs
    ):
        coords = {
            "temperature": xarray.Variable(
                "temperature",
                temperature.value,
                attrs={"unit": temperature.unit.to_string()},
            ),
            "wavelength": xarray.Variable(
                "wavelength",
                wavelength.value,
                attrs={"unit": wavelength.unit.to_string()},
            ),
        }
        if density is None:
            coords["density"] = xarray.Variable(
                "temperature", density.value, attrs={"unit": density.unit.to_string()}
            )
        self._da = xarray.DataArray(
            spectra.data,
            dims=["temperature", "wavelength"],
            coords=coords,
            attrs={"unit": spectra.unit.to_string(), **kwargs},
        )

    def __repr__(self):
        return self._da.__repr__()

    @property
    @u.quantity_input
    def temperature(self) -> u.K:
        return u.Quantity(self._da.temperature.data, self._da.temperature.attrs["unit"])

    @property
    @u.quantity_input
    def wavelength(self) -> u.Angstrom:
        return u.Quantity(self._da.wavelength.data, self._da.wavelength.attrs["unit"])

    @property
    @u.quantity_input
    def data(self) -> u.photon * u.cm**3 / (u.s * u.AA * u.sr):
        return u.Quantity(self._da.data, self._da.attrs["unit"])


class AbstractChannel(abc.ABC):
    """
    An abstract base class for defining instrument channels.

    .. caution::

       This abstract class is still under development and may change
       in the near future.
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
        weighted by the wavelength response of the instrument. For more information,
        see the `topic guide on instrument response <>`_.

        Parameters
        ----------
        source_spectra: `SourceSpectra`
        obstime: any parsed by `~sunpy.time.parse_time`, optional
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
        return (
            self.geometrical_area
            * self.mirror_reflectance
            * self.filter_transmittance
            * self.quantum_efficiency
            * self.degradation(obstime=obstime)
        )

    @property
    @u.quantity_input
    def energy_per_photon(self) -> u.eV / u.photon:
        return self.wavelength.to("eV", equivalencies=u.spectral()) / u.photon

    @abc.abstractmethod
    @u.quantity_input
    def degradation(self, obstime=None):
        ...

    @abc.abstractproperty
    @u.quantity_input
    def geometrical_area(self) -> u.cm**2:
        ...

    @abc.abstractproperty
    @u.quantity_input
    def mirror_reflectance(self) -> u.dimensionless_unscaled:
        ...

    @abc.abstractproperty
    @u.quantity_input
    def filter_transmittance(self) -> u.dimensionless_unscaled:
        ...

    @abc.abstractproperty
    @u.quantity_input
    def quantum_efficiency(self) -> u.dimensionless_unscaled:
        ...

    @abc.abstractproperty
    @u.quantity_input
    def camera_gain(self) -> u.DN / u.electron:
        ...

    @abc.abstractproperty
    @u.quantity_input
    def energy_per_electron(self) -> u.eV / u.electron:
        ...

    @abc.abstractproperty
    @u.quantity_input
    def pixel_solid_angle(self) -> u.steradian / u.pixel:
        ...

    @abc.abstractproperty
    @u.quantity_input
    def wavelength(self) -> u.Angstrom:
        ...
