import numpy as np
import pytest

import astropy.units as u

from sunkit_instruments.response import SourceSpectra, get_temperature_response
from sunkit_instruments.response.abstractions import AbstractChannel, EmissionModel


class TestChannel(AbstractChannel):
    @property
    @u.quantity_input
    def wavelength(self) -> u.Angstrom:
        return np.linspace(100, 200, 100) * u.AA

    @u.quantity_input
    def degradation(self, obstime=None) -> u.dimensionless_unscaled:
        return 1.0

    @property
    @u.quantity_input
    def geometrical_area(self) -> u.cm**2:
        return 10 * u.cm**2

    @property
    @u.quantity_input
    def mirror_reflectance(self) -> u.dimensionless_unscaled:
        return np.exp(
            -((self.wavelength - self.wavelength[0]) / self.wavelength[0]).decompose()
        )

    @property
    @u.quantity_input
    def filter_transmittance(self) -> u.dimensionless_unscaled:
        return np.exp(
            -((self.wavelength - 150 * u.AA) ** 2 / (1 * u.AA) ** 2).decompose()
        )

    @property
    @u.quantity_input
    def effective_quantum_efficiency(self) -> u.dimensionless_unscaled:
        return np.ones(self.wavelength.shape)

    @property
    @u.quantity_input
    def camera_gain(self) -> u.DN / u.electron:
        return 2 * u.DN / u.electron

    @property
    @u.quantity_input
    def energy_per_electron(self) -> u.eV / u.electron:
        return 3.65 * u.eV / u.electron

    @property
    @u.quantity_input
    def pixel_solid_angle(self) -> u.steradian / u.pixel:
        return (1 * u.arcsec) ** 2 / u.pixel


@pytest.fixture
def fake_channel():
    return TestChannel()


def test_effective_area(fake_channel):
    assert isinstance(fake_channel.effective_area(), u.Quantity)


def test_wavelength_response(fake_channel):
    assert isinstance(fake_channel.wavelength_response(), u.Quantity)


@pytest.fixture
def fake_spectra():
    temperature = np.logspace(5, 8, 100) * u.K
    density = 1e15 * u.cm ** (-3) * u.K / temperature
    wavelength = np.linspace(50, 250, 1000) * u.AA
    data = np.random.rand(*temperature.shape + wavelength.shape) * u.Unit(
        "photon cm3 s-1 sr-1 Angstrom-1"
    )
    return SourceSpectra(temperature, wavelength, data, density=density)


def test_spectra_repr(fake_spectra):
    assert isinstance(fake_spectra.__repr__(), str)


@pytest.mark.parametrize('obstime', [None, '2020-01-01'])
def test_temperature_response(fake_channel, fake_spectra, obstime):
    temp_response = fake_spectra.temperature_response(fake_channel, obstime=obstime)
    assert isinstance(temp_response, u.Quantity)
    assert temp_response.shape == fake_spectra.temperature.shape


class TimeDependentChannel(TestChannel):
    @u.quantity_input
    def degradation(self, obstime=None) -> u.dimensionless_unscaled:
        return 2.0 if obstime is not None else 1.0


class FakeEmissionModel:
    temperature = np.array([1, 2]) * u.MK

    def calculate_temperature_response(self, channel):
        return np.repeat(channel.wavelength_response().sum(), self.temperature.size)


def test_fake_emission_model_conforms_to_protocol():
    assert isinstance(FakeEmissionModel(), EmissionModel)
    assert not isinstance(object(), EmissionModel)


def test_get_temperature_response_with_emission_model():
    temperature, temp_response = get_temperature_response(
        TimeDependentChannel(), FakeEmissionModel()
    )
    assert isinstance(temperature, u.Quantity)
    assert isinstance(temp_response, u.Quantity)
    assert temp_response.shape == temperature.shape


def test_get_temperature_response_with_emission_model_obstime():
    channel = TimeDependentChannel()
    model = FakeEmissionModel()
    _, response_without_obstime = get_temperature_response(channel, model)
    _, response_with_obstime = get_temperature_response(
        channel, model, obstime="2020-01-01"
    )
    assert np.allclose(response_with_obstime.value, 2 * response_without_obstime.value)


def test_get_temperature_response_with_unsupported_source(fake_channel):
    with pytest.raises(TypeError, match="temperature_response"):
        get_temperature_response(fake_channel, object())
