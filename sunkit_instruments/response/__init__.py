"""
A subpackage for computing instrument responses
"""

from sunkit_instruments.response.abstractions import AbstractChannel, EmissionModel, LineEmissionModel
from sunkit_instruments.response.linelist import (
    LineListEmissionModel,
    chianti_line_list,
    get_line_list,
    line_list_cache_path,
    line_list_from_emission_model,
)
from sunkit_instruments.response.spectral import create_response_function, get_spectral_response
from sunkit_instruments.response.thermal import SourceSpectra, get_temperature_response
