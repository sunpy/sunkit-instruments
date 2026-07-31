"""
========================================================
GOES flare detection comparison between 1s, 1min and HEK
========================================================

This example shows how to use `sunkit_instruments` to detect flares in
GOES-XRS data, compares the results between the high-cadence (1-second)
and 1-minute averaged data product, and compare both against the
officially reported GOES flare list from the HEK.

We will look at an active periodn in May 2024 with many large flare inxlucig
an X-class and also a quiet period during solar minimum, in May 2020, where
the largest flare in the windowis only A-class.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import astropy.units as u
from astropy.time import Time

from sunpy import timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import TimeRange

from sunkit_instruments.goes_xrs import (
    find_flares_naive,
    find_goes_flares,
    flareclass_to_flux,
    get_goes_event_list,
)

###############################################################################
# We start by searching for the GOES XRS data for our chosen time range. NOAA
# provides two relevant data products: ``flx1s``, the high-cadence (1-second
# for GOES-R satellites) flux, and ``avg1m``, the flux already averaged to a
# 1-minute cadence.
#
# NOAA's own flare lists are built from a primary satellite, falling back to a
# secondary satellite whenever the primary's data is flagged as poor quality
# or missing, so we fetch both satellites here and do the same. Since we
# repeat this fetch-and-merge step for two different time ranges (and pairs
# of satellites) below, we wrap it in a function.


def merge_primary_secondary(primary_df, secondary_df):
    primary_good = primary_df[primary_df["xrsb_quality"] == 0]
    all_minutes = pd.Index(primary_df.index.floor("min").unique())
    good_minutes = pd.Index(primary_good.index.floor("min").unique())
    bad_minutes = all_minutes.difference(good_minutes)

    secondary_good = secondary_df[secondary_df["xrsb_quality"] == 0]
    fill = secondary_good[secondary_good.index.floor("min").isin(bad_minutes)]

    return pd.concat([primary_good["xrsb"], fill["xrsb"]]).sort_index()


def fetch_merged_fluxes(tr, satellites):
    query = Fido.search(
        tr,
        a.Instrument.xrs,
        a.goes.SatelliteNumber(satellites[0]) | a.goes.SatelliteNumber(satellites[1]),
        a.Resolution.flx1s | a.Resolution.avg1m,
    )
    files = Fido.fetch(query)

    dataframes = {}
    for satellite in satellites:
        for resolution, key in [("avg1m", "1min"), ("flx1s", "1s")]:
            matches = [f for f in files if f"g{satellite}_" in f and resolution in f]
            dataframes[(satellite, key)] = ts.TimeSeries(matches, concatenate=True).to_dataframe()

    merged_1min = merge_primary_secondary(
        dataframes[(satellites[0], "1min")], dataframes[(satellites[1], "1min")]
    )
    merged_1s = merge_primary_secondary(
        dataframes[(satellites[0], "1s")], dataframes[(satellites[1], "1s")]
    )
    return merged_1s, merged_1min


###############################################################################
# `~sunkit_instruments.goes_xrs.find_goes_flares` implements the NOAA SWPC
# algorithm for detecting flares directly from a `~sunpy.timeseries.sources.XRSTimeSeries`,
# and classifies each detection into a GOES class. It always rebins the input
# onto a 1-minute cadence internally (matching the definition of the
# algorithm), regardless of the native cadence of the data, so we expect both
# timeseries to agree closely.


def to_flares(series, min_class, **kwargs):
    df = series.to_frame(name="xrsb")
    goes_ts = ts.TimeSeries(df, {}, {"xrsb": u.W / u.m**2}, source="xrs")
    return find_goes_flares(goes_ts, min_class=min_class, **kwargs)


###############################################################################
# To compare detections against the official list, we stack three panels,
# one each for the 1-second detections, the 1-minute detections, and the
# official HEK list. Each panel shows the same long channel flux (the merged
# 1-minute data, for all three, so the panels are directly comparable) with
# every flare shaded from its start to its end time, and its start, peak and
# end times marked with different line styles. This is wrapped in a function
# so we can reuse it below.

vline_styles = {"start_time": ":", "peak_time": "-", "end_time": "--"}


def plot_comparison(merged_1min, flares_1s, flares_1min, hek_events, title):
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 9))

    panels = [
        (axes[0], flares_1s, "1-second data", "tab:blue"),
        (axes[1], flares_1min, "1-minute averaged data", "tab:orange"),
        (axes[2], hek_events, "Official GOES event list (HEK)", "tab:red"),
    ]

    for ax, events, label, color in panels:
        ax.plot(merged_1min.index, merged_1min.to_numpy(), color="black")
        for event in events:
            ax.axvspan(
                event["start_time"].datetime, event["end_time"].datetime, color=color, alpha=0.2
            )
            for name, style in vline_styles.items():
                ax.axvline(event[name].datetime, color=color, linestyle=style, alpha=0.8)
        ax.set_yscale("log")
        ax.set_ylabel(label)

    for name, style in vline_styles.items():
        axes[0].plot([], [], color="black", linestyle=style, label=name.replace("_", " "))
    axes[0].legend(loc="upper right")

    axes[-1].set_xlabel("Time")
    fig.suptitle(title)
    fig.autofmt_xdate()
    plt.show()


###############################################################################
# Finally, we look at exactly which HEK events our algorithm failed to detect.
# For each one, we plot the long channel flux in a window around its reported
# peak time, together with its start/peak/end times, and shade in any of our
# own detected events (from the 1-minute data) that overlap the same window.
#
# We also try to diagnose *why* each one was missed, by walking back from the
# HEK peak along the monotonic rise leading up to it: if that rise is shorter
# than ``min_rise_time`` it is annotated "too short", and if it doesn't reach
# ``min_rise_ratio`` times its own starting flux it is annotated "not steep".
# This is again wrapped in a function so we can reuse it below.

MIN_RISE_TIME = 4  # minutes; matches find_flares' default min_rise_time
MIN_RISE_RATIO = 1.4  # matches find_flares' default min_rise_ratio


def diagnose_miss(series, peak_time):
    flux = series.to_numpy()
    idx = series.index.get_indexer([peak_time], method="nearest")[0]
    j = idx
    while j > 0 and flux[j - 1] < flux[j]:
        j -= 1
    duration = (series.index[idx] - series.index[j]) / pd.Timedelta(minutes=1)
    ratio = flux[idx] / flux[j] if flux[j] > 0 else np.inf
    if duration < MIN_RISE_TIME:
        return "too short"
    if ratio < MIN_RISE_RATIO:
        return "too shollow"
    return None


def plot_missed_events(hek_events, flares_1min, merged_1min, window=pd.Timedelta(minutes=25)):
    our_peaks = flares_1min["peak_time"].datetime64
    missed_events = [
        event
        for event in hek_events
        if np.abs((our_peaks - event["peak_time"].datetime64) / np.timedelta64(1, "m")).min() > 5
    ]
    print(f"{len(missed_events)} of {len(hek_events)} HEK events were not detected by our algorithm")

    if not missed_events:
        return missed_events

    ncols = 4
    nrows = int(np.ceil(len(missed_events) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.4 * nrows), squeeze=False)

    for ax, event in zip(axes.flat, missed_events):
        peak = pd.Timestamp(event["peak_time"].datetime)
        lo, hi = peak - window, peak + window
        segment = merged_1min.loc[lo:hi]
        ax.plot(segment.index, segment.to_numpy(), color="black", lw=1)
        ax.set_xlim(lo, hi)

        for row in flares_1min:
            if row["start_time"].datetime <= hi and row["end_time"].datetime >= lo:
                ax.axvspan(
                    max(row["start_time"].datetime, lo),
                    min(row["end_time"].datetime, hi),
                    color="tab:blue",
                    alpha=0.15,
                )
        for name, style in vline_styles.items():
            ax.axvline(event[name].datetime, color="tab:red", linestyle=style, alpha=0.9)

        nearest_idx = int(
            np.argmin(np.abs((our_peaks - event["peak_time"].datetime64) / np.timedelta64(1, "m")))
        )
        ax.axvline(
            flares_1min["peak_time"][nearest_idx].datetime,
            color="steelblue",
            linestyle=vline_styles["peak_time"],
            alpha=0.9,
        )

        ax.set_yscale("log")
        title = f"{event['goes_class']} {peak:%m-%d %H:%M}"
        diagnosis = diagnose_miss(merged_1min, peak)
        if diagnosis is not None:
            title += f"\n({diagnosis})"
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", labelsize=6, rotation=30)
        ax.tick_params(axis="y", labelsize=7)

    for ax in axes.flat[len(missed_events):]:
        ax.axis("off")

    fig.suptitle(
        "HEK flares not detected by our algorithm\n"
        "(red = HEK start/peak/end; steelblue = our nearest peak; "
        "blue shading = a nearby detection of ours, if any)"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()
    return missed_events


###############################################################################
# The other direction matters too: events we detect that HEK doesn't report
# at all. These aren't necessarily wrong, HEK's automated list has gaps of
# its own, but a cluster of them is worth a look, especially in a quiet
# period where a much lower detection threshold makes us far more sensitive
# to noise and instrumental artifacts (see the electron contamination
# correction discussed above).


def plot_extra_events(hek_events, flares_1min, merged_1min, window=pd.Timedelta(minutes=25)):
    hek_peaks = np.array([event["peak_time"].datetime64 for event in hek_events])
    extra_events = [
        row
        for row in flares_1min
        if len(hek_peaks) == 0
        or np.abs((hek_peaks - row["peak_time"].datetime64) / np.timedelta64(1, "m")).min() > 5
    ]
    print(f"{len(extra_events)} of {len(flares_1min)} detections have no matching HEK event")

    if not extra_events:
        return extra_events

    ncols = 4
    nrows = int(np.ceil(len(extra_events) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.4 * nrows), squeeze=False)

    for ax, row in zip(axes.flat, extra_events):
        peak = pd.Timestamp(row["peak_time"].datetime)
        lo, hi = peak - window, peak + window
        segment = merged_1min.loc[lo:hi]
        ax.plot(segment.index, segment.to_numpy(), color="black", lw=1)
        ax.set_xlim(lo, hi)

        for event in hek_events:
            if event["start_time"].datetime <= hi and event["end_time"].datetime >= lo:
                ax.axvspan(
                    max(event["start_time"].datetime, lo),
                    min(event["end_time"].datetime, hi),
                    color="tab:red",
                    alpha=0.15,
                )
        for name, style in vline_styles.items():
            ax.axvline(row[name].datetime, color="steelblue", linestyle=style, alpha=0.9)

        ax.set_yscale("log")
        ax.set_title(f"{row['goes_class']} {peak:%m-%d %H:%M}", fontsize=9)
        ax.tick_params(axis="x", labelsize=6, rotation=30)
        ax.tick_params(axis="y", labelsize=7)

    for ax in axes.flat[len(extra_events):]:
        ax.axis("off")

    fig.suptitle(
        "Our detections with no matching HEK event\n"
        "(steelblue = our start/peak/end; red shading = a nearby HEK event, if any)"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()
    return extra_events


###############################################################################
# `~sunkit_instruments.goes_xrs.find_goes_flares` reverse engineers several
# refinements on top of NOAA's own 3-rule description of start, peak and end
# (an extending-window search for the rise ratio, a peak-locking mechanism,
# and noise-tolerant patching of an interrupted decay's end time) to better
# match the official list. `~sunkit_instruments.goes_xrs.find_flares_naive`
# implements *only* those 3 rules exactly as described, in a single pass,
# with none of that extra machinery, so we can compare the two against each
# other below.


def match_fraction(flares, hek_events, tol_min=5):
    peaks = flares["peak_time"].datetime64
    matched = sum(
        1
        for event in hek_events
        if np.abs((peaks - event["peak_time"].datetime64) / np.timedelta64(1, "m")).min() <= tol_min
    )
    return matched, len(hek_events)


###############################################################################
# Active period
# -------------
# GOES-16 was the primary satellite throughout May 2024, with GOES-18 as its
# secondary.

tr = a.Time("2024-05-10 00:00", "2024-05-15 23:59")
# if offline run for entire moth or longer
# tr = a.Time("2024-05-01 00:00", "2024-05-31 23:59")

merged_1s, merged_1min = fetch_merged_fluxes(tr, satellites=(16, 18))

MIN_CLASS = "B1.0"

hek_events = get_goes_event_list(TimeRange(tr.start, tr.end), goes_class_filter=MIN_CLASS)

print("Flares reported in the official GOES event list (from the HEK):", len(hek_events))

###############################################################################
# Quiet period
# ------------
# Now let's set up a quiet period too, during solar minimum. May 2020 was
# otherwise extremely quiet, but NOAA AR 12765/12766, one of the first
# flaring active regions of Solar Cycle 25, produced a short burst of
# activity (including an M1.1, the largest flare of the month) between
# 2020-05-27 and 2020-05-29, so we focus there. GOES-16 was still the primary
# satellite in May 2020, but its secondary at the time was GOES-17 rather
# than GOES-18. Since we're now looking for much weaker events than the
# active period above, we also lower ``min_class`` and ``flux_threshold``
# (the quiet-Sun background floor used to detect a rise at all) from the
# B1.0 defaults down to A1.0.
#
# (some of the extra events will dispear if use the extended info in the 1min avg files)

tr_quiet = a.Time("2020-05-26 00:00", "2020-05-31 23:59")
# if offline run for entire month or longer
# tr_quiet = a.Time("2020-05-01 00:00", "2020-05-31 23:59")

merged_1s_quiet, merged_1min_quiet = fetch_merged_fluxes(tr_quiet, satellites=(16, 17))

MIN_CLASS_QUIET = "A1.0"

hek_events_quiet = get_goes_event_list(
    TimeRange(tr_quiet.start, tr_quiet.end), goes_class_filter=MIN_CLASS_QUIET
)

print("Flares reported in the official GOES event list (from the HEK):", len(hek_events_quiet))

###############################################################################
# Naive method
# ------------
# Let's first see how the naive, 3-rule-only method does against the
# official list, on both periods, using the 1-minute averaged flux for each
# (the algorithm is defined for 1-minute averages).

naive_flares = find_flares_naive(
    Time(merged_1min.index.values, format="datetime64"),
    merged_1min.to_numpy() * u.W / u.m**2,
    min_flux=flareclass_to_flux(MIN_CLASS),
    flux_threshold=flareclass_to_flux(MIN_CLASS),
)
naive_matched, total = match_fraction(naive_flares, hek_events)
print(
    f"Naive (3-rule) method, active period: {len(naive_flares)} flares detected, "
    f"{naive_matched}/{total} HEK events matched"
)

naive_flares_quiet = find_flares_naive(
    Time(merged_1min_quiet.index.values, format="datetime64"),
    merged_1min_quiet.to_numpy() * u.W / u.m**2,
    min_flux=flareclass_to_flux(MIN_CLASS_QUIET),
    flux_threshold=flareclass_to_flux(MIN_CLASS_QUIET),
)
naive_matched_quiet, total_quiet = match_fraction(naive_flares_quiet, hek_events_quiet)
print(
    f"Naive (3-rule) method, quiet period:  {len(naive_flares_quiet)} flares detected, "
    f"{naive_matched_quiet}/{total_quiet} HEK events matched"
)

###############################################################################
# Improved method
# ----------------
# Now let's see how `~sunkit_instruments.goes_xrs.find_goes_flares` does on
# the same two periods. Both cadences typically agree closely with each
# other, and recover substantially more of the events in the official list
# than the naive method above. It only implements the *automated* part of
# NOAA's detection algorithm though, and most of the remaining gap comes from
# local sub-peaks that occur during a single, still-elevated, complex flare:
# NOAA's list reports these as their own separate events, while our peak,
# once locked in, only starts a new flare after flux has clearly turned over
# for a sustained period. A smaller share of the gap is made up of slow,
# gradual events whose rise is never steep enough to trip the automated
# trigger at all, which the official list can include via manual entries
# from SWPC operators.
#
# To see this, we stack three panels, one each for the 1-second detections,
# the 1-minute detections, and the official HEK list, and separately plot
# the HEK events our algorithm missed.

flares_1s = to_flares(merged_1s, min_class=MIN_CLASS)
flares_1min = to_flares(merged_1min, min_class=MIN_CLASS)

print("Flares detected from 1-second data:", len(flares_1s))
print("Flares detected from 1-minute averaged data:", len(flares_1min))

refined_matched, _ = match_fraction(flares_1min, hek_events)
print(
    f"Improved method, active period:       {len(flares_1min)} flares detected, "
    f"{refined_matched}/{total} HEK events matched"
)

###############################################################################
# Now let see a comparison

plot_comparison(
    merged_1min,
    flares_1s,
    flares_1min,
    hek_events,
    "Long channel flux [W/m$^2$] with detected/reported flares",
)


###############################################################################
# Missed events compare to HEK

missed_events = plot_missed_events(hek_events, flares_1min, merged_1min)

###############################################################################
# Extra events we detect that HEK doesn't report

extra_events = plot_extra_events(hek_events, flares_1min, merged_1min)

flares_1s_quiet = to_flares(
    merged_1s_quiet, min_class=MIN_CLASS_QUIET, flux_threshold=flareclass_to_flux(MIN_CLASS_QUIET)
)
flares_1min_quiet = to_flares(
    merged_1min_quiet, min_class=MIN_CLASS_QUIET, flux_threshold=flareclass_to_flux(MIN_CLASS_QUIET)
)

print("Flares detected from 1-second data:", len(flares_1s_quiet))
print("Flares detected from 1-minute averaged data:", len(flares_1min_quiet))

refined_matched_quiet, _ = match_fraction(flares_1min_quiet, hek_events_quiet)
print(
    f"Improved method, quiet period:        {len(flares_1min_quiet)} flares detected, "
    f"{refined_matched_quiet}/{total_quiet} HEK events matched"
)

###############################################################################
# Now let see a comparison

plot_comparison(
    merged_1min_quiet,
    flares_1s_quiet,
    flares_1min_quiet,
    hek_events_quiet,
    "Long channel flux [W/m$^2$] with detected/reported flares (quiet period)",
)

###############################################################################
# Missed events compare to HEK

missed_events_quiet = plot_missed_events(hek_events_quiet, flares_1min_quiet, merged_1min_quiet)

###############################################################################
# Extra events we detect that HEK doesn't report

extra_events_quiet = plot_extra_events(hek_events_quiet, flares_1min_quiet, merged_1min_quiet)
