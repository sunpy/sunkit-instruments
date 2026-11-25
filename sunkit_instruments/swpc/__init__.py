"""
This module provides a Fido client and parser for the NOAA SWPC "Edited Solar
Events List" files.
"""
from sunkit_instruments.swpc.client import SWPCEventsClient
from sunkit_instruments.swpc.parser import parse_swpc_event_archive, parse_swpc_event_file, parse_swpc_events

__all__ = ["SWPCEventsClient", "parse_swpc_events", "parse_swpc_event_file", "parse_swpc_event_archive"]
