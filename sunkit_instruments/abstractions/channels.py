"""This module defines abstractions related to instrument channels."""
import abc

import numpy as np
import xarray
from scipy.interpolate import interp1d

import astropy.units as u

__all__ = ["AbstractChannel", "SourceSpectra"]


class SourceSpectra:
    def __init__(self, temperature, wavelength, spectra, density=None, **kwargs):
        # Internal xarray representation
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

    @property
    def temperature(self) -> u.K:
        return u.Quantity(self._da.temperature.data, self._da.temperature.attrs["unit"])

    @property
    def wavelength(self) -> u.Angstrom:
        return u.Quantity(self._da.wavelength.data, self._da.wavelength.attrs["unit"])

    @property
    def data(self) -> u.photon * u.cm**3 / (u.s * u.AA * u.sr):
        return u.Quantity(self._da.data, self._da.attrs["unit"])


class AbstractChannel(abc.ABC):
    """An abstract base class for defining instrument channels.

    .. caution::

       This abstract class is still under development and may change
       in the near future.
    """

    def temperature_response(
        self, source_spectra, obstime=None, **kwargs
    ) -> u.cm**5 * u.DN / (u.pix * u.s):
        wave_response = self.wavelength_response(obstime=obstime, **kwargs)
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
            source_spectra * wave_response_interp * np.gradient(self.wavelength)
        ).sum(axis=1)
        return source_spectra.temperature, temperature_response

    def wavelength_response(
        self, obstime=None, **kwargs
    ) -> u.cm**2 * u.DN * u.steradian / (u.photon * u.pixel):
        area_eff = self.effective_area(obstime=obstime, **kwargs)
        return (
            area_eff
            * self.photon_to_energy
            * self.pixel_area
            / (self.energy_to_electron * self.camera_gain)
        )

    def effective_area(self, obstime=None) -> u.cm**2:
        return (
            self.geometrical_area
            * self.mirror_reflectance
            * self.filter_transmittance
            * self.quantum_efficiency
            * self.degradation(obstime=obstime)
        )

    @property
    def photon_to_energy(self) -> u.eV / u.photon:
        return self.wavelength.to("eV", equivalencies=u.spectral()) / u.photon

    @abc.abstractmethod
    def degradation(self, obstime=None):
        ...

    @abc.abstractproperty
    def geometrical_area(self) -> u.cm**2:
        ...

    @abc.abstractproperty
    def mirror_reflectance(self):
        ...

    @abc.abstractproperty
    def filter_transmittance(self):
        ...

    @abc.abstractproperty
    def quantum_efficiency(self):
        ...

    @abc.abstractproperty
    def camera_gain(self) -> u.electron / u.DN:
        ...

    @abc.abstractproperty
    def energy_to_electron(self) -> u.electron / u.eV:
        ...

    @abc.abstractproperty
    def pixel_area(self) -> u.cm**2 / u.pixel:
        ...

    @abc.abstractproperty
    def wavelength(self) -> u.Angstrom:
        ...
