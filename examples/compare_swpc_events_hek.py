"""
=====================================================
Comparing the SWPC events list to the HEK GOES flares
=====================================================

This example shows how to search for, download and parse a NOAA SWPC
"Edited Solar Events List" using `~sunkit_instruments.swpc`, and compares the
GOES X-ray flares (``XRA``) it contains to the equivalent flare list held by
the `Heliophysics Event Knowledgebase (HEK) <https://www.lmsal.com/hek/>`__,
alongside the underlying GOES X-ray flux the flares were detected in.
"""

import bisect

import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table

from sunpy import timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import TimeRange

import sunkit_instruments.swpc  # noqa: F401
from sunkit_instruments.goes_xrs import flareclass_to_flux
from sunkit_instruments.swpc import parse_swpc_events

###############################################################################
# We define a helper function which, for a given time range, searches for and
# downloads the SWPC events list using
# `sunkit_instruments.swpc.SWPCEventsClient`, parses out the GOES X-ray
# events (event type ``"XRA"``) with `~sunkit_instruments.swpc.parse_swpc_events`,
# and gets the equivalent flare list straight from the HEK.
# Flares in the HEK reported by the SWPC (i.e. with a feature recognition
# method, or "FRM", of ``"SWPC"``) are themselves derived from this same
# events list.
# Simply importing `sunkit_instruments.swpc` is enough to register the client
# with `~sunpy.net.Fido`.
#
# It then plots the GOES class of each flare (converted to X-ray flux with
# `~sunkit_instruments.goes_xrs.flareclass_to_flux`) from both catalogs
# against the underlying GOES X-ray flux, itself downloaded and read in with
# `sunpy.timeseries.TimeSeries`.
#
# Depending on how old ``timerange`` is, each search result may be a single
# day's file or a whole year's archive (see
# `~sunkit_instruments.swpc.SWPCEventsClient`); `~sunkit_instruments.swpc.parse_swpc_events`
# handles both transparently, and the ``time_range`` argument keeps a yearly
# archive from contributing more than what was actually asked for.
# Each file is parsed separately, and only the flares with a real (unmasked)
# ``max_time`` are kept, rather than combining the parsed tables with
# `~astropy.table.vstack`: ``vstack`` does not preserve the mask on a
# `~astropy.time.Time` column, so doing this would silently turn any missing
# ``max_time`` into a bogus, but seemingly valid, date.


def plot_swpc_and_hek_flares(timerange, satellite_number=16):
    swpc_results = Fido.search(timerange, a.Instrument.swpc_events)
    swpc_files = Fido.fetch(swpc_results[0])
    time_range = TimeRange(timerange.start, timerange.end)
    swpc_flare_rows = []
    for file in swpc_files:
        table = parse_swpc_events(file, time_range=time_range)
        xra = table[table["event_type"] == "XRA"]
        xra["goes_class"] = [row["particulars"].split()[0] for row in xra]
        swpc_flare_rows.extend(row for row in xra if not row["max_time"].mask)
    swpc_flare_rows.sort(key=lambda row: row["max_time"])
    # Built with ``Table(rows=...)`` rather than `~astropy.table.vstack`, since
    # ``vstack`` does not preserve the mask on a `~astropy.time.Time` column.
    swpc_flares = Table(rows=swpc_flare_rows, names=swpc_flare_rows[0].colnames)

    hek_results = Fido.search(timerange, a.hek.EventType("FL"), a.hek.FRM.Name == "SWPC")
    hek_flares = hek_results["hek"]

    goes_results = Fido.search(timerange, a.Instrument.xrs, a.Resolution.avg1m, a.goes.SatelliteNumber(satellite_number))
    goes_ts = ts.TimeSeries(Fido.fetch(goes_results), concatenate=True)

    # A HEK flare is considered "matched" if its peak time falls within 2
    # minutes of a SWPC event's peak time.
    has_match = np.array([
        any(abs((flare["event_peaktime"] - row["max_time"]).sec) < 120 for row in swpc_flares)
        for flare in hek_flares
    ])
    hek_matched, hek_unmatched = hek_flares[has_match], hek_flares[~has_match]

    print(f"{len(swpc_flares)} XRA events in the SWPC events list.")
    print(f"{len(hek_flares)} flares in the HEK, {len(hek_unmatched)} with no matching SWPC event.")
    for flare in hek_unmatched:
        peak = flare["event_peaktime"]
        idx = bisect.bisect_left(swpc_flares["max_time"], peak)
        before, after = max(idx - 1, 0), min(idx + 2, len(swpc_flares))
        print(
            f"Unmatched HEK flare: {flare['event_peaktime'].iso} peak, {flare['fl_goescls']} "
            "(would fall between the SWPC events below):"
        )
        print(swpc_flares[before:after], '\n')

    fig, ax = plt.subplots()
    goes_ts.to_dataframe()["xrsb"].plot(
        ax=ax, color="grey", lw=0.75, label=f"GOES-{satellite_number} XRS (1-8 Å)",
    )
    for i, row in enumerate(swpc_flares):
        ax.scatter(
            row["max_time"].datetime.item(), flareclass_to_flux(row["goes_class"]).value,
            marker="x", color="C0", s=80, label="SWPC events list" if i == 0 else None,
        )
    for i, flare in enumerate(hek_matched):
        ax.scatter(
            flare["event_peaktime"].datetime, flareclass_to_flux(flare["fl_goescls"]).value,
            marker="+", color="C1", s=100, label="HEK (SWPC)" if i == 0 else None,
        )
    for i, flare in enumerate(hek_unmatched):
        ax.scatter(
            flare["event_peaktime"].datetime, flareclass_to_flux(flare["fl_goescls"]).value,
            marker="s", color="C3", s=50, label="HEK (SWPC), no SWPC match" if i == 0 else None,
        )
    ax.set_yscale("log")
    ax.set_ylabel("GOES X-ray flux [W/m$^2$]")
    ax.set_xlabel("Time [UTC]")
    ax.legend()
    fig.autofmt_xdate()
    return fig


###############################################################################
# A busy period
# -------------
# First, let's look at 2024 May 10-15, the week of the X3.9 flare (from NOAA
# active region 13664) that triggered the severe "Mother's Day" geomagnetic
# storm. The two catalogs agree closely, including on the X3.9 flare itself,
# although the HEK list contains 3 flares with no matching entry in the SWPC
# edited list. Two of these are short brightenings whose end time coincides
# exactly with the start of the following, larger flare, suggesting they were
# folded into that event when the "duty forecaster" manually edited the
# list (see the SWPC events list `README
# <ftp://ftp.swpc.noaa.gov/pub/indices/events/README>`__). The third is an
# X3.4 flare that is clearly present in the underlying GOES flux, but has no
# corresponding entry at all in the SWPC edited list - a genuine miss, of the
# kind the README itself acknowledges are possible.

plot_swpc_and_hek_flares(a.Time("2024-05-10", "2024-05-15 23:59"))

###############################################################################
# A quiet period
# ---------------
# Now, let's compare the same catalogs over 2020 May, a much quieter month
# deep in solar minimum. Nearly every day has no reported events at all,
# aside from a single active day (2020 May 29) which both catalogs agree on
# exactly, down to a single, small M1.1 flare.

plot_swpc_and_hek_flares(a.Time("2020-05-01", "2020-05-31 23:59"))

plt.show()
