import numpy as np
import pandas as pd

import astropy.units as u
from astropy.table import QTable
from astropy.time import Time

from sunpy import timeseries as ts

from .goes_xrs import flareclass_to_flux, flux_to_flareclass

__all__ = ["find_flares_naive", "find_flares", "find_goes_flares"]


def find_flares_naive(
    times,
    flux,
    min_flux=1e-6 * u.W / u.m**2,
    flux_threshold=1e-7 * u.W / u.m**2,
    min_rise_time=4 * u.min,
    min_rise_ratio=1.4,
    min_decay_time=4 * u.min,
    decay_fraction=0.5,
    cadence=1 * u.min,
):
    """
    Detect flares in a timeseries using a literal reading of NOAA's 3-rule
    description, with none of the extra refinements `find_flares` adds to
    better match the official GOES flare list:

    * A flare starts at the first sample, ``i``, where flux is at least
      ``flux_threshold`` and rises, uninterrupted, for exactly
      ``min_rise_time``, reaching a value more than ``min_rise_ratio`` times
      the flux at ``i`` within that fixed window. Unlike `find_flares`, the
      ratio is not given any further samples to build towards if it is not
      met within this window alone.
    * The flare peaks at the maximum flux reached before the flux falls for
      ``min_decay_time`` consecutive samples, exactly as in `find_flares`.
    * The flare ends once the flux has decayed ``decay_fraction`` of the way
      back down from the peak flux to the flux at the start of the rise. If
      that level is never reached before the data ends, or before a later,
      unrelated flare's own decay happens to cross it, the end time is
      whatever sample the plain scan first satisfies this at, rather than
      being patched against the next flare's start.

    This function exists to demonstrate the effect of those refinements: it
    recovers substantially fewer of the officially reported GOES flares than
    `find_flares` does, and is not intended for production use.

    The algorithm is defined for 1-minute averaged flux, so ``times`` and
    ``flux`` are first rebinned onto the requested ``cadence`` (by averaging)
    regardless of their native sampling.

    Parameters
    ----------
    times : `~astropy.time.Time`
        The time of each measurement.
    flux : `~astropy.units.Quantity`
        The measured signal to search for flares in, in any unit. NaN values
        are treated as data gaps and are never considered part of a rise or
        decay.
    min_flux : `~astropy.units.Quantity`, optional
        The minimum peak flux for a flare to be included in the returned table.
        Must be convertible to the unit of ``flux``. Defaults to 1e-6 W/m**2
        (i.e. GOES class C1). Pass `None` to disable filtering.
    flux_threshold : `~astropy.units.Quantity`, optional
        The minimum flux for a rise to be considered above the quiet-Sun
        background at all. Must be convertible to the unit of ``flux``.
        Defaults to 1e-7 W/m**2 (i.e. GOES class B1).
    min_rise_time : `~astropy.units.Quantity`, optional
        The exact duration over which the rise ratio must be met. Defaults to
        4 minutes. Internally converted to a number of rebinned samples using
        ``cadence``, so must be at least two ``cadence``.
    min_rise_ratio : `float`, optional
        The minimum ratio between the flux at the end and at the start of the
        rise window for it to be considered a genuine flare. Defaults to 1.4.
    min_decay_time : `~astropy.units.Quantity`, optional
        The duration of consecutively decreasing flux that marks the end of
        the rise to peak flux. Defaults to 4 minutes. Internally converted to
        a number of rebinned samples using ``cadence``, so must be at least
        one ``cadence``.
    decay_fraction : `float`, optional
        The fraction of the way back down from the peak flux to the pre-flare
        background flux that defines the end of a flare. Defaults to 0.5.
    cadence : `~astropy.units.Quantity`, optional
        The cadence that ``times`` and ``flux`` are rebinned onto (by
        averaging) before the algorithm is run. Defaults to 1 minute.

    Returns
    -------
    `~astropy.table.QTable`
        A table of the detected flares with columns ``start_time``, ``peak_time``,
        ``end_time`` and ``peak_flux``.

    See Also
    --------
    find_flares : the refined version of this algorithm used elsewhere in
        this package, which recovers substantially more of the officially
        reported GOES flares.

    References
    ----------
    * `GOES XRS Guide <https://ngdc.noaa.gov/stp/satellite/goes/doc/GOES_XRS_readme.pdf>`__
    """
    if len(times) != len(flux):
        raise ValueError(
            f"times and flux must be the same length, got {len(times)} and {len(flux)}."
        )

    if not isinstance(flux, u.Quantity):
        raise TypeError(f"flux must be an astropy Quantity, not {type(flux)}")
    flux_unit = flux.unit
    flux = flux.value

    if min_flux is not None and not isinstance(min_flux, u.Quantity):
        raise TypeError(f"min_flux must be an astropy Quantity or None, not {type(min_flux)}")
    if not isinstance(flux_threshold, u.Quantity):
        raise TypeError(f"flux_threshold must be an astropy Quantity, not {type(flux_threshold)}")

    threshold_flux = -np.inf if min_flux is None else min_flux.to_value(flux_unit)
    detection_floor = flux_threshold.to_value(flux_unit)

    n_rise_samples = round((min_rise_time / cadence).to_value(u.dimensionless_unscaled))
    n_decay_samples = round((min_decay_time / cadence).to_value(u.dimensionless_unscaled))
    if n_rise_samples < 2 or n_decay_samples < 1:
        raise ValueError(
            "min_rise_time must be at least two cadences, and min_decay_time at least one."
        )

    series = pd.Series(flux, index=pd.DatetimeIndex(times.datetime64))
    series = series.resample(pd.Timedelta(cadence.to_value(u.s), unit="s")).mean()

    times = Time(series.index.values, format="datetime64")
    flux = series.to_numpy()
    n = len(flux)

    start_indices, peak_indices, end_indices, peak_fluxes = [], [], [], []

    i = 0
    while i <= n - n_rise_samples:
        window = flux[i : i + n_rise_samples]
        if np.any(np.isnan(window)) or window[0] < detection_floor:
            i += 1
            continue
        if not np.all(np.diff(window) > 0) or not (window[-1] > min_rise_ratio * window[0]):
            i += 1
            continue

        # `i` is the pre-flare background level. As in `find_flares`, the
        # peak is locked in as soon as flux has been falling for
        # `n_decay_samples` consecutive samples.
        background = flux[i]
        peak_index = i
        peak_flux = flux[i]
        consecutive_decreases = 0
        previous_flux = peak_flux

        j = i + 1
        while j < n:
            current_flux = flux[j]
            if np.isnan(current_flux):
                consecutive_decreases = 0
                j += 1
                continue

            if current_flux > peak_flux:
                peak_flux = current_flux
                peak_index = j
                consecutive_decreases = 0
            elif current_flux < previous_flux:
                consecutive_decreases += 1
                if consecutive_decreases >= n_decay_samples:
                    break
            else:
                consecutive_decreases = 0

            previous_flux = current_flux
            j += 1

        # Unlike `find_flares`, a plain scan for the first sample that
        # reaches `decay_level`: no check for whether a later uptick is a
        # genuine new rise, and no patching against the next flare's start if
        # that level is never reached.
        decay_level = background + (peak_flux - background) * (1 - decay_fraction)
        end_index = n - 1
        for k in range(peak_index + 1, n):
            if not np.isnan(flux[k]) and flux[k] <= decay_level:
                end_index = k
                break

        if peak_flux >= threshold_flux:
            start_indices.append(i)
            peak_indices.append(peak_index)
            end_indices.append(end_index)
            peak_fluxes.append(peak_flux)

        i = peak_index + 1

    return QTable(
        {
            "start_time": times[start_indices],
            "peak_time": times[peak_indices],
            "end_time": times[end_indices],
            "peak_flux": u.Quantity(peak_fluxes, flux_unit),
        }
    )


def _rise_starts_at(flux, i, n, detection_floor, min_rise_ratio, n_rise_samples):
    """
    Whether a genuine rise starts at index `i`: an uninterrupted monotonic
    increase, staying above `detection_floor`, that is at least
    `n_rise_samples` long and at some point exceeds `min_rise_ratio` times
    the flux at `i`.
    """
    if np.isnan(flux[i]) or flux[i] < detection_floor:
        return False
    j = i + 1
    while j < n:
        current_flux = flux[j]
        if np.isnan(current_flux) or current_flux <= flux[j - 1]:
            return False
        if j - i + 1 >= n_rise_samples and current_flux > min_rise_ratio * flux[i]:
            return True
        j += 1
    return False


def find_flares(
    times,
    flux,
    min_flux=1e-6 * u.W / u.m**2,
    flux_threshold=1e-7 * u.W / u.m**2,
    min_rise_time=4 * u.min,
    min_rise_ratio=1.4,
    min_decay_time=4 * u.min,
    decay_fraction=0.5,
    cadence=1 * u.min,
):
    """
    Detect flares in a timeseries.

    This implements a reverse engineered version of NOAA SWPC algorithm used
    to build the official GOES flare lists:

    * A flare starts at the first sample, ``i``, of an uninterrupted run of
      increasing flux above ``flux_threshold`` that is at least
      ``min_rise_time`` long and, at some point along the way, reaches a
      value more than ``min_rise_ratio`` times the flux at ``i``. Real flares
      often build gradually, so this ratio need not be met within the first
      ``min_rise_time`` alone; as long as the rise stays uninterrupted, later
      samples still count towards the same candidate start.
    * The flare peaks at the maximum flux reached before the flux falls for
      ``min_decay_time`` consecutive samples: once that many consecutive
      samples of decrease are seen, the peak is locked in, and a later,
      stronger re-brightening starts its own, separate flare rather than
      being folded into this one.
    * The flare ends once the flux has decayed ``decay_fraction`` of the way
      back down from the peak flux to the flux at the start of the rise (0.5,
      i.e. half way back to the pre-flare background, by default). If a new
      rise begins before that level is reached, the decay is considered
      interrupted, and this flare's end is set to that next flare's start
      instead, matching the convention used in the official GOES event lists.

    The algorithm is defined for 1-minute averaged flux, so ``times`` and
    ``flux`` are first rebinned onto the requested ``cadence`` (by averaging)
    regardless of their native sampling.

    Notes
    -----
    This is the automated detection algorithm; the official NOAA/GOES event
    list also includes events added manually by SWPC operators for edge cases
    the automated algorithm cannot handle. This function will not find those.

    Parameters
    ----------
    times : `~astropy.time.Time`
        The time of each measurement.
    flux : `~astropy.units.Quantity`
        The measured signal to search for flares in, in any unit. This need
        not be a physically calibrated flux; a detector count rate or data
        number works just as well, as long as ``min_flux`` and
        ``flux_threshold`` are given in a compatible unit. NaN values are
        treated as data gaps and are never considered part of a rise or decay.
    min_flux : `~astropy.units.Quantity`, optional
        The minimum peak flux for a flare to be included in the returned table.
        This only filters the returned table; it does not affect detection
        (see ``flux_threshold`` for that). Must be convertible to the unit of
        ``flux``. Defaults to 1e-6 W/m**2 (i.e. GOES class C1). Pass `None`
        to disable filtering.
    flux_threshold : `~astropy.units.Quantity`, optional
        The minimum flux for samples to be considered above the quiet-Sun
        background at all. Must be convertible to the unit of ``flux``.
        Defaults to 1e-7 W/m**2 (i.e. GOES class B1), matching the NOAA
        algorithm.
    min_rise_time : `~astropy.units.Quantity`, optional
        The duration of consecutively increasing flux that marks the start of a
        flare. Defaults to 4 minutes. Internally converted to a number of
        rebinned samples using ``cadence``, so must be at least two ``cadence``.
    min_rise_ratio : `float`, optional
        The minimum ratio between the flux at the end and at the start of the
        rise window for it to be considered a genuine flare, rather than a
        marginal fluctuation in the background. Defaults to 1.4, matching the
        NOAA algorithm.
    min_decay_time : `~astropy.units.Quantity`, optional
        The duration of consecutively decreasing flux that marks the end of the
        rise to peak flux. Defaults to 4 minutes. Internally converted to a
        number of rebinned samples using ``cadence``, so must be at least one
        ``cadence``.
    decay_fraction : `float`, optional
        The fraction of the way back down from the peak flux to the pre-flare
        background flux that defines the end of a flare. Defaults to 0.5,
        i.e. the flare ends when the flux decays half way back to the
        background level, matching the NOAA algorithm.
    cadence : `~astropy.units.Quantity`, optional
        The cadence that ``times`` and ``flux`` are rebinned onto (by averaging)
        before the algorithm is run. Defaults to 1 minute, matching the standard
        NOAA algorithm.

    Returns
    -------
    `~astropy.table.QTable`
        A table of the detected flares with columns ``start_time``, ``peak_time``,
        ``end_time`` and ``peak_flux``.

    References
    ----------
    * `GOES XRS Guide <https://ngdc.noaa.gov/stp/satellite/goes/doc/GOES_XRS_readme.pdf>`__
    """
    if len(times) != len(flux):
        raise ValueError(
            f"times and flux must be the same length, got {len(times)} and {len(flux)}."
        )

    if not isinstance(flux, u.Quantity):
        raise TypeError(f"flux must be an astropy Quantity, not {type(flux)}")
    flux_unit = flux.unit
    flux = flux.value

    if min_flux is not None and not isinstance(min_flux, u.Quantity):
        raise TypeError(f"min_flux must be an astropy Quantity or None, not {type(min_flux)}")
    if not isinstance(flux_threshold, u.Quantity):
        raise TypeError(f"flux_threshold must be an astropy Quantity, not {type(flux_threshold)}")

    threshold_flux = -np.inf if min_flux is None else min_flux.to_value(flux_unit)
    detection_floor = flux_threshold.to_value(flux_unit)

    n_rise_samples = round((min_rise_time / cadence).to_value(u.dimensionless_unscaled))
    n_decay_samples = round((min_decay_time / cadence).to_value(u.dimensionless_unscaled))
    if n_rise_samples < 2 or n_decay_samples < 1:
        raise ValueError(
            "min_rise_time must be at least two cadences, and min_decay_time at least one."
        )

    series = pd.Series(flux, index=pd.DatetimeIndex(times.datetime64))
    series = series.resample(pd.Timedelta(cadence.to_value(u.s), unit="s")).mean()

    times = Time(series.index.values, format="datetime64")
    flux = series.to_numpy()
    n = len(flux)

    # Collect every candidate event first, regardless of `min_flux`, so that
    # an interrupted flare's end can always be patched against the very next
    # flare found (see below), not just the next one that happens to pass
    # the `min_flux` filter.
    all_start_indices, all_peak_indices, all_end_indices = [], [], []
    all_peak_fluxes, all_end_confirmed = [], []

    i = 0
    while i <= n - n_rise_samples:
        # Real flares often build gradually, so `min_rise_ratio` is frequently
        # not met within the first `n_rise_samples` samples alone; as long as
        # the rise from `i` is uninterrupted, later samples still count
        # towards the same candidate start at `i`.
        if not _rise_starts_at(flux, i, n, detection_floor, min_rise_ratio, n_rise_samples):
            i += 1
            continue

        # `i` is the pre-flare background level. The peak is locked in as soon
        # as flux has been falling for `n_decay_samples` consecutive samples:
        # a later, stronger re-brightening starts its own, separate flare
        # rather than being folded into this one.
        background = flux[i]
        peak_index = i
        peak_flux = flux[i]
        consecutive_decreases = 0
        previous_flux = peak_flux

        j = i + 1
        while j < n:
            current_flux = flux[j]
            if np.isnan(current_flux):
                consecutive_decreases = 0
                j += 1
                continue

            if current_flux > peak_flux:
                peak_flux = current_flux
                peak_index = j
                consecutive_decreases = 0
            elif current_flux < previous_flux:
                consecutive_decreases += 1
                if consecutive_decreases >= n_decay_samples:
                    break
            else:
                consecutive_decreases = 0

            previous_flux = current_flux
            j += 1

        # The flare ends once flux has decayed `decay_fraction` of the way
        # back down from the (now locked) peak to the pre-flare background.
        # A local uptick alone doesn't interrupt this search (that would make
        # the end far too sensitive to ordinary noise in the decay); it only
        # does so if it is the start of a genuine new rise, in which case the
        # end is filled in below, once that next flare's own start is known.
        decay_level = background + (peak_flux - background) * (1 - decay_fraction)
        end_index = None
        previous_flux = peak_flux
        previous_index = peak_index
        k = peak_index + 1
        while k < n:
            current_flux = flux[k]
            if np.isnan(current_flux):
                k += 1
                continue
            if current_flux <= decay_level:
                end_index = k
                break
            if current_flux > previous_flux and _rise_starts_at(
                flux, previous_index, n, detection_floor, min_rise_ratio, n_rise_samples
            ):
                break
            previous_flux = current_flux
            previous_index = k
            k += 1

        all_start_indices.append(i)
        all_peak_indices.append(peak_index)
        all_end_indices.append(end_index if end_index is not None else n - 1)
        all_end_confirmed.append(end_index is not None)
        all_peak_fluxes.append(peak_flux)

        i = peak_index + 1

    # A flare whose decay was interrupted by the next one starting shares its
    # end time with that next flare's start, rather than an arbitrary or
    # missing value; this also matches the convention used in the official
    # GOES event lists.
    for idx in range(len(all_end_indices) - 1):
        if not all_end_confirmed[idx]:
            all_end_indices[idx] = all_start_indices[idx + 1]

    start_indices, peak_indices, end_indices, peak_fluxes = [], [], [], []
    for start, peak, end, peak_flux in zip(
        all_start_indices, all_peak_indices, all_end_indices, all_peak_fluxes
    ):
        if peak_flux >= threshold_flux:
            start_indices.append(start)
            peak_indices.append(peak)
            end_indices.append(end)
            peak_fluxes.append(peak_flux)

    return QTable(
        {
            "start_time": times[start_indices],
            "peak_time": times[peak_indices],
            "end_time": times[end_indices],
            "peak_flux": u.Quantity(peak_fluxes, flux_unit),
        }
    )


def find_goes_flares(goes_ts, min_class="C1", **kwargs):
    """
    Detect flares in a GOES XRS timeseries.

    This extracts the long channel (``xrsb``) flux of ``goes_ts`` and runs it
    through `~sunkit_instruments.goes_xrs.find_flares`, which implements a version
    of the NOAA SWPC algorithm used to build the official GOES flare lists,
    then classifies each detected flare's peak flux into a GOES class.

    Parameters
    ----------
    goes_ts : `~sunpy.timeseries.sources.XRSTimeSeries`
        The GOES XRS timeseries containing the long channel (``xrsb``) flux (in W/m**2).
    min_class : `str`, optional
        The minimum GOES class for a flare to be included in the returned table,
        e.g. "C1", "M1". Defaults to "C1".
    **kwargs :
        Additional keyword arguments are passed to `~sunkit_instruments.goes_xrs.find_flares`.

    Returns
    -------
    `~astropy.table.QTable`
        A table of the detected flares with columns ``start_time``, ``peak_time``,
        ``end_time``, ``goes_class`` and ``peak_flux``.

    References
    ----------
    * `GOES XRS Guide <https://ngdc.noaa.gov/stp/satellite/goes/doc/GOES_XRS_readme.pdf>`__
    """
    if not isinstance(goes_ts, ts.XRSTimeSeries):
        raise TypeError(
            f"Input time series must be a XRSTimeSeries instance, not {type(goes_ts)}"
        )
    if "xrsb" not in goes_ts.columns:
        raise ValueError("The input time series does not contain a 'xrsb' column.")

    df = goes_ts.to_dataframe()
    longflux = df["xrsb"].copy()
    if "xrsb_quality" in df.columns:
        longflux[df["xrsb_quality"] != 0] = np.nan

    times = Time(longflux.index.values, format="datetime64")
    flux = u.Quantity(longflux.to_numpy(), "W/m**2")

    if "min_flux" in kwargs:
        raise TypeError("min_flux is derived from min_class; pass min_class instead.")

    flares = find_flares(times, flux, min_flux=flareclass_to_flux(min_class), **kwargs)
    flares["goes_class"] = [flux_to_flareclass(peak_flux) for peak_flux in flares["peak_flux"]]
    return flares["start_time", "peak_time", "end_time", "goes_class", "peak_flux"]
