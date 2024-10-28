"""
Classes for computing the temperature response
"""
import xarray

import astropy.units as u

__all__ = ["SourceSpectra", "get_temperature_response"]


def get_temperature_response(channel, spectra, obstime=None):
    """
    Calculate the temperature response function for a given instrument channel
    and input spectra.

    Parameters
    ----------
    channel: `~sunkit_instruments.response.abstractions.AbstractChannel`
    spectra: `~sunkit_instruments.response.SourceSpectra`
    obstime: any format parsed by `sunpy.time.parse_time`, optional

    Returns
    -------
    temperature: `~astropy.units.Quantity`
    response: `~astropy.units.Quantity`

    See Also
    --------
    sunkit_instruments.response.SourceSpectra.temperature_response
    """
    return spectra.temperature, spectra.temperature_response(channel, obstime=obstime)


class SourceSpectra:
    """
    Source spectra as a function of temperature and wavelength.

    The source spectra describes how a plasma emits (under the optically-thin
    assumption) as a function of both temperature and wavelength. This source
    spectra is typically computed using a database like CHIANTI by summing the
    emission spectra of many ions as well as the continuum emission. For more
    information, see the topic guide on instrument response.

    Parameters
    ----------
    temperature: `~astropy.units.Quantity`
        1D array describing the variation along the temperature axis.
    wavelength: `~astropy.units.Quantity`
        1D array describing the variation along the wavelength axis.
    spectra: `~astropy.units.Quantity`
        Source spectra as a 2D array. The first axis should correspond to temperature and the
        second axis should correspond to wavelength.
    density: `~astropy.units.Quantity`, optional
        1D array describing the variation in density along the temperature axis. It is assumed
        that temperature and density are dependent.
    meta: `dict`, optional
        Any optional metadata to attach to the spectra, e.g. abundance model, CHIANTI version.
    """

    @u.quantity_input
    def __init__(
        self,
        temperature: u.K,
        wavelength: u.Angstrom,
        spectra: u.photon * u.cm**3 / (u.s * u.Angstrom * u.steradian),
        density: u.cm ** (-3) = None,
        meta=None,
    ):
        self.meta = meta
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
            attrs={"unit": spectra.unit.to_string(), **self.meta},
        )

    def __repr__(self):
        return self._da.__repr__()

    def __str__(self):
        return self._da.__str__()

    def _repr_html_(self):
        return self._da._repr_html_()

    @property
    def meta(self):
        return self._meta

    @meta.setter
    def meta(self, x):
        if x is None:
            self._meta = {}
        elif isinstance(x, dict):
            self._mata = x
        else:
            raise TypeError(f'Unsupported metadata type {type(x)}')

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
    def data(self) -> u.photon * u.cm**3 / (u.s * u.Angstrom * u.steradian):
        return u.Quantity(self._da.data, self._da.attrs["unit"])

    @u.quantity_input
    def temperature_response(
        self, channel, obstime=None
    ) -> u.cm**5 * u.DN / (u.pixel * u.s):
        """
        Temperature response function for a given instrument channel.

        The temperature response function describes the sensitivity of an imaging
        instrument as a function of temperature. The temperature response is
        calculated by integrating the source spectra over the wavelength dimension,
        weighted by the wavelength response of the instrument.

        Parameters
        ----------
        channel: `~sunkit_instruments.response.abstractions.AbstractChannel`
            The relevant instrument channel object used to compute the wavelength
            response function.
        obstime: any format parsed by `sunpy.time.parse_time`, optional
            A time of a particular observation. This is used to calculated any
            time-dependent instrument degradation.
        """
        wave_response = channel.wavelength_response(obstime=obstime)
        da_response = xarray.DataArray(
            wave_response.to_value(wave_response.unit),
            dims=['wavelength'],
            coords=[xarray.Variable('wavelength',
                                   channel.wavelength.to_value('Angstrom'),
                                   attrs={'unit': 'Angstrom'}),],
            attrs={'unit': wave_response.unit.to_string()}
        )
        spec_interp = self._da.interp(
            wavelength=da_response.wavelength,
            kwargs={'bounds_error': False, 'fill_value': 0.0},
        )
        spec_interp_weighted = spec_interp * da_response
        temp_response = spec_interp_weighted.integrate(coord='wavelength')
        final_unit = (u.Unit(spec_interp.unit)
                      * u.Unit(da_response.attrs['unit'])
                      * u.Unit(spec_interp_weighted.wavelength.unit))
        return u.Quantity(temp_response.data, final_unit)
