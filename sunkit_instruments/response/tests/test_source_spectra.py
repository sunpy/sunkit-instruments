"""
Smoke tests for SourceSpectra
"""
import numpy as np
import pytest

import astropy.units as u

from sunkit_instruments.response import SourceSpectra


@pytest.mark.parametrize(('density', 'meta'), [
    (None, None),
    (1e15*u.K/u.cm**3, None),
    (None, {'abundance_model': 'test'}),
    (1e15*u.K/u.cm**3, {'abundance_model': 'test'}),
])
def test_create_source_spectra(density, meta):
    temperature = np.logspace(4, 9, 100) * u.K
    wavelength = np.linspace(1, 1000, 50) * u.Angstrom
    if density is not None:
        density = density / temperature
    data_shape = temperature.shape + wavelength.shape
    data = np.random.rand(*data_shape) * u.Unit("photon cm3 s-1 sr-1 Angstrom-1")
    spec = SourceSpectra(
        temperature,
        wavelength,
        data,
        density=density,
        meta=meta,
    )
    assert isinstance(spec.meta, dict)
    assert isinstance(spec.data, u.Quantity)
    assert isinstance(spec.wavelength, u.Quantity)
    assert isinstance(spec.temperature, u.Quantity)
    if density is not None:
        assert isinstance(spec.density, u.Quantity)
    else:
        with pytest.raises(ValueError, match="No density data available."):
            spec.density
