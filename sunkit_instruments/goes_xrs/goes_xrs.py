import numpy as np

from astropy import units as u

from sunpy.time import parse_time

GOES_CONVERSION_DICT = {
    "X": u.Quantity(1e-4, "W/m^2"),
    "M": u.Quantity(1e-5, "W/m^2"),
    "C": u.Quantity(1e-6, "W/m^2"),
    "B": u.Quantity(1e-7, "W/m^2"),
    "A": u.Quantity(1e-8, "W/m^2"),
}

__all__ = [
    "get_goes_event_list",
    "flux_to_flareclass",
    "flareclass_to_flux",
]


def get_goes_event_list(timerange, goes_class_filter=None):
    """
    Retrieve list of flares detected by GOES within a given time range.

    Parameters
    ----------
    timerange : `sunpy.time.TimeRange`
        The time range to download the event list for.
    goes_class_filter: `str`, optional
        A string specifying a minimum GOES class for inclusion in the list,
        e.g., "M1", "X2".

    Returns
    -------
    `list`:
        A list of all the flares found for the given time range.
    """
    # Importing hek here to avoid calling code that relies on optional dependencies.
    from sunpy.net import attrs, hek

    # use HEK module to search for GOES events
    client = hek.HEKClient()
    event_type = "FL"
    tstart = timerange.start
    tend = timerange.end

    # query the HEK for a list of events detected by the GOES instrument
    # between tstart and tend (using a GOES-class filter)
    if goes_class_filter:
        result = client.search(
            attrs.Time(tstart, tend),
            attrs.hek.EventType(event_type),
            attrs.hek.FL.GOESCls > goes_class_filter,
            attrs.hek.OBS.Observatory == "GOES",
        )
    else:
        result = client.search(
            attrs.Time(tstart, tend),
            attrs.hek.EventType(event_type),
            attrs.hek.OBS.Observatory == "GOES",
        )

    # want to condense the results of the query into a more manageable
    # dictionary
    # keep event data, start time, peak time, end time, GOES-class,
    # location, active region source (as per GOES list standard)
    # make this into a list of dictionaries
    goes_event_list = []

    for r in result:
        goes_event = {
            "event_date": parse_time(r["event_starttime"]).strftime("%Y-%m-%d"),
            "start_time": parse_time(r["event_starttime"]),
            "peak_time": parse_time(r["event_peaktime"]),
            "end_time": parse_time(r["event_endtime"]),
            "goes_class": str(r["fl_goescls"]),
            "goes_location": (r["event_coord1"], r["event_coord2"]),
            "noaa_active_region": r["ar_noaanum"],
        }
        goes_event_list.append(goes_event)

    return goes_event_list


def flareclass_to_flux(flareclass):
    """
    Converts a GOES flare class into the corresponding X-ray flux.

    Parameters
    ----------
    flareclass : str
        The case-insensitive flare class (e.g., 'X3.2', 'm1.5', 'A9.6').

    Returns
    -------
    flux : `~astropy.units.Quantity`
        X-ray flux between 1 and 8 Angstroms as measured near Earth in W/m^2.

    Raises
    ------
    TypeError
        Input must be a string.

    Examples
    --------
    >>> from sunkit_instruments.goes_xrs import flareclass_to_flux
    >>> flareclass_to_flux('A1.0')
    <Quantity 1.e-08 W / m2>
    >>> flareclass_to_flux('c4.7')
    <Quantity 4.7e-06 W / m2>
    >>> flareclass_to_flux('X2.4')
    <Quantity 0.00024 W / m2>
    """
    if not isinstance(flareclass, str):
        raise TypeError(f"Input must be a string, not {type(flareclass)}")
    # TODO should probably make sure the string is in the expected format.
    flareclass = flareclass.upper()
    # invert the conversion dictionary
    # conversion_dict = {v: k for k, v in GOES_CONVERSION_DICT.items()}
    return float(flareclass[1:]) * GOES_CONVERSION_DICT[flareclass[0]]


@u.quantity_input
def flux_to_flareclass(goesflux: u.watt / u.m**2):
    """
    Converts X-ray flux into the corresponding GOES flare class.

    Parameters
    ----------
    flux : `~astropy.units.Quantity`
        X-ray flux between 1 and 8 Angstroms (usually measured by GOES) as
        measured at the Earth in W/m^2

    Returns
    -------
    flareclass : str
        The flare class e.g.: 'X3.2', 'M1.5', 'A9.6'.

    Raises
    ------
    ValueError
        Flux cannot be negative.

    References
    ----------
    `Solar Flare Classification <https://en.wikipedia.org/wiki/Solar_flare#Classification>`_

    Examples
    --------
    >>> from sunkit_instruments.goes_xrs import flux_to_flareclass
    >>> import astropy.units as u
    >>> flux_to_flareclass(1e-08 * u.watt/u.m**2)
    'A1'
    >>> flux_to_flareclass(4.7e-06 * u.watt/u.m**2)
    'C4.7'
    >>> flux_to_flareclass(0.00024 * u.watt/u.m**2)
    'X2.4'
    >>> flux_to_flareclass(7.8e-09 * u.watt/u.m**2)
    'A0.78'
    >>> flux_to_flareclass(0.00682 * u.watt/u.m**2)
    'X68.2'
    """

    if goesflux.value < 0:
        raise ValueError("Flux cannot be negative")

    decade = np.floor(np.log10(goesflux.to("W/m**2").value))
    # invert the conversion dictionary
    conversion_dict = {v: k for k, v in GOES_CONVERSION_DICT.items()}
    if decade < -8:
        str_class = "A"
        decade = -8
    elif decade > -4:
        str_class = "X"
        decade = -4
    else:
        str_class = conversion_dict.get(u.Quantity(10**decade, "W/m**2"))
    goes_subclass = 10**-decade * goesflux.to("W/m**2").value
    return f"{str_class}{goes_subclass:.3g}"
