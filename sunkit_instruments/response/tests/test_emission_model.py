"""
Tests for the emission-model line-list adapter and get_spectral_response,
using a fake model that mimics the synthesizAR EmissionModel surface.
No fiasco or synthesizAR required.
"""

import numpy as np
import pytest
import xarray as xr

import astropy.constants as const
import astropy.units as u

from sunkit_instruments.response import (
    LineListEmissionModel,
    get_spectral_response,
    get_temperature_response,
    line_list_from_emission_model,
)
from sunkit_instruments.response.abstractions import EmissionModel, LineEmissionModel
from sunkit_instruments.response.tests.test_channel import TestChannel

TEMPERATURE = np.array([5e5, 1e6, 2e6]) * u.K
DENSITY = np.array([1e9]) * u.cm**-3


class FakeIon:
    def __init__(self, ion_name, ion_name_roman, atomic_number, ionization_fraction):
        self.ion_name = ion_name
        self.ion_name_roman = ion_name_roman
        self.atomic_number = atomic_number
        self.ionization_fraction = ionization_fraction


class FakeLineEmissionModel:
    temperature = TEMPERATURE
    density = DENSITY

    def __init__(self):
        self.fe9 = FakeIon("fe_9", "Fe IX", 26, np.array([0.1, 0.8, 0.2]))
        self.missing = FakeIon("na_3", "Na III", 11, np.ones(3))
        self._ions = [self.fe9, self.missing]
        self.emissivity = u.Quantity(np.ones((3, 1, 2)), "cm3 ph s-1")
        self.wavelengths = np.array([174.5, 171.073]) * u.AA

    def __iter__(self):
        return iter(self._ions)

    def get_line_emissivity(self, ion, transition=None):
        if ion is self.missing:
            # placeholder for ions without level-population data
            return [0] * u.AA, u.Quantity(np.zeros((3, 1, 1)), "cm3 ph s-1")
        return self.wavelengths, self.emissivity


def test_fake_model_conforms_to_protocol():
    assert isinstance(FakeLineEmissionModel(), LineEmissionModel)


def test_adapter_schema_and_sorting():
    line_list = line_list_from_emission_model(FakeLineEmissionModel())
    assert line_list.sizes["trans_index"] == 2  # missing-data ion skipped
    np.testing.assert_allclose(line_list.wavelength.values, [171.073, 174.5])  # sorted
    assert list(line_list.full_name.values) == ["Fe IX 171.073", "Fe IX 174.500"]
    assert list(line_list.atomic_number.values) == [26, 26]
    assert line_list.gofnt.dims == ("logT", "density", "abundance", "trans_index")
    np.testing.assert_allclose(line_list.logT.values, np.log10(TEMPERATURE.to_value(u.K)))
    assert u.Unit(line_list.gofnt.attrs["units"]) == u.Unit("erg cm3 / (s sr)")


def test_adapter_gofnt_values():
    """gofnt = emissivity x ionization fraction x photon energy / 4 pi sr."""
    model = FakeLineEmissionModel()
    line_list = line_list_from_emission_model(model)
    photon_energy = (const.h * const.c / (171.073 * u.AA)).to_value(u.erg)
    expected = 1.0 * model.fe9.ionization_fraction * photon_energy / (4 * np.pi)
    gofnt_171 = line_list.gofnt.sel(trans_index=line_list.wavelength == 171.073).squeeze()
    np.testing.assert_allclose(gofnt_171.values, expected, rtol=1e-12)


def test_adapter_can_defer_ionization_fraction():
    model = FakeLineEmissionModel()
    with_ioneq = line_list_from_emission_model(model)
    without_ioneq = line_list_from_emission_model(model, include_ionization_fraction=False)
    ratio = (with_ioneq.gofnt / without_ioneq.gofnt).isel(density=0, abundance=0, trans_index=0)
    np.testing.assert_allclose(ratio.values, model.fe9.ionization_fraction)


def test_adapter_empty_model_raises():
    model = FakeLineEmissionModel()
    model._ions = [model.missing]
    with pytest.raises(ValueError, match="no line emissivities"):
        line_list_from_emission_model(model)


def test_get_spectral_response_from_model():
    channel = TestChannel()
    response = get_spectral_response(
        channel,
        FakeLineEmissionModel(),
        vdop=np.array([-100.0, 0.0, 100.0]) * u.km / u.s,
        wavelength_range=[170.0, 176.0],
        num_lines_keep=1,
    )
    assert {"line", "logT", "vdop", "density", "abundance", "wavelength"} <= set(response.response.dims)
    lines = [str(line) for line in response.line.values]
    assert lines[0] == "Fe IX 171.073"
    assert any("remaining" in line.lower() for line in lines[1:])
    assert "erg" in response.response.attrs["units"]


def test_get_spectral_response_accepts_line_list_dataset():
    channel = TestChannel()
    line_list = line_list_from_emission_model(FakeLineEmissionModel())
    response = get_spectral_response(
        channel,
        line_list,
        wavelength_range=[170.0, 176.0],
        num_lines_keep=2,
    )
    assert response.sizes["line"] == 2


@pytest.fixture
def chianti_style_line_list():
    logT = np.array([5.8, 6.0, 6.2])
    gofnt = np.array([1.0, 5.0, 2.0]) * 1e-25
    line_list = xr.Dataset(
        {
            "wavelength": ("trans_index", np.array([171.073])),
            "atomic_number": ("trans_index", np.array([26])),
            "gofnt": (
                ("logT", "pressure", "abundance", "trans_index"),
                gofnt[:, np.newaxis, np.newaxis, np.newaxis],
            ),
            "full_name": ("trans_index", np.array(["Fe IX 171.073"], dtype=object)),
        },
        coords={"logT": logT, "pressure": [3e15], "abundance": ["fake"]},
    )
    line_list.gofnt.attrs["units"] = "erg cm3 / (s sr)"
    return line_list


def test_line_list_emission_model_conforms_to_protocol(chianti_style_line_list):
    assert isinstance(LineListEmissionModel(chianti_style_line_list), EmissionModel)


def test_line_list_emission_model_temperature_response(chianti_style_line_list):
    channel = TestChannel()
    model = LineListEmissionModel(chianti_style_line_list)
    temperature, response = get_temperature_response(channel, model)

    assert u.allclose(temperature, 10 ** chianti_style_line_list.logT.data * u.K)
    assert response.unit.physical_type == u.Unit("cm5 DN / (pix s)").physical_type
    # single line at 171.073: K(T) = gofnt(T) * R(171.073) * photons/erg
    wave_response = channel.wavelength_response()
    r_at_line = np.interp(171.073, channel.wavelength.to_value(u.AA), wave_response.value)
    photons_per_erg = 1 / (const.h * const.c / (171.073 * u.AA)).to_value(u.erg)
    expected = chianti_style_line_list.gofnt.values.squeeze() * r_at_line * photons_per_erg
    np.testing.assert_allclose(response.to_value("cm5 DN / (pix s)"), expected, rtol=1e-12)


def test_line_list_emission_model_zero_outside_channel(chianti_style_line_list):
    line_list = chianti_style_line_list.copy()
    line_list["wavelength"] = ("trans_index", np.array([500.0]))  # outside 100-200 A channel
    _, response = get_temperature_response(TestChannel(), LineListEmissionModel(line_list))
    np.testing.assert_array_equal(response.value, 0.0)


def test_line_list_emission_model_requires_single_grid():
    line_list = xr.Dataset(
        {
            "wavelength": ("trans_index", np.array([171.073])),
            "gofnt": (("logT", "pressure", "trans_index"), np.ones((2, 2, 1))),
        },
        coords={"logT": [5.8, 6.0], "pressure": [3e15, 3e16]},
    )
    with pytest.raises(ValueError, match="single pressure"):
        LineListEmissionModel(line_list)
