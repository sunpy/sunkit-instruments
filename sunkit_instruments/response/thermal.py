"""
Classes for computing the temperature response
"""
import xarray

import astropy.units as u

__all__ = ["SourceSpectra", "TemperatureResponse"]


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
    temperature: `~astropy.units.Quantity`
        1D array describing the variation along the temperature axis
    wavelength: `~astropy.units.Quantity`
        1D array describing the variation along the wavelength axis
    spectra: `~astropy.units.Quantity`
        Source spectra as a 2D array. The first axis should correspond to temperature and the second axis should correspond to wavelength.
    density: `~astropy.units.Quantity`, optional
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
        if density is not None:
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


class TemperatureResponse:
    ...
