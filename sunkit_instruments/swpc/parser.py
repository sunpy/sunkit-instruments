"""
This module provides a parser for the NOAA SWPC "Edited Solar Events List" files.

See the `README
<ftp://ftp.swpc.noaa.gov/pub/indices/events/README>`__ for a full description
of the file format.
"""
import re
import tarfile
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np

from astropy.table import Column, MaskedColumn, QTable, vstack
from astropy.time import Time

__all__ = ["parse_swpc_events", "parse_swpc_event_file", "parse_swpc_event_archive"]

# Column byte offsets within the fixed 80 character wide data rows.
# These were derived from, and validated against, the column headers and
# real files spanning from 1997 to the present day.
_COLUMNS = {
    "event_id": (0, 4),
    "flag": (4, 7),
    "begin": (10, 15),
    "max": (17, 22),
    "end": (27, 32),
    "observatory": (33, 37),
    "quality": (38, 40),
    "event_type": (41, 46),
    "location_or_frequency": (47, 57),
    "particulars": (57, 74),
    "region_number": (74, 80),
}

# An arbitrary placeholder for missing `max_time`/`end_time` values, chosen
# to be unambiguously distinct from any real SWPC report date (which starts
# in 1996) while still being new enough that ERFA does not consider it a
# "dubious year".
_MISSING_TIME_SENTINEL = Time("1980-01-01T00:00:00")

_NO_EVENTS = re.compile(r"^\s*NO EVENT REPORTS\.?\s*$", re.IGNORECASE)
_DATE_LINE = re.compile(r"^:Date:\s*(\d{4})\s+(\d{1,2})\s+(\d{1,2})\s*$")
_EDITED_EVENTS_LINE = re.compile(
    r"Edited Events for\s+(\d{4})\s+([A-Za-z]{3})\s+(\d{1,2})", re.IGNORECASE
)


def _get_report_date(lines):
    """
    Find the UTC day that a SWPC events list file reports on.

    This is deliberately not the ``:Created:`` timestamp (or, in the older
    file format, the timestamp on the second line of the file) as that is
    when the file was written out, which can be a number of days after the
    UTC day that the file's events belong to.
    """
    for line in lines:
        match = _DATE_LINE.match(line)
        if match:
            year, month, day = (int(value) for value in match.groups())
            return datetime(year, month, day)
        match = _EDITED_EVENTS_LINE.search(line)
        if match:
            year, month_abbr, day = match.groups()
            return datetime.strptime(f"{year} {month_abbr} {day}", "%Y %b %d")
    raise ValueError("Could not find the report date in the SWPC events file.")


def _parse_time(raw, report_date, *, after=None):
    """
    Parse a ``HHMM`` field into a full `~astropy.time.Time`.

    Parameters
    ----------
    raw : `str`
        The raw field, e.g. ``"0318"``, ``"A0146"`` or ``"////"``.
    report_date : `datetime.datetime`
        The UTC day that the file (and therefore the ``Begin`` time) belongs to.
    after : `datetime.datetime`, optional
        If the parsed time is earlier than this, roll over to the next UTC day.
        This is used for the ``Max`` and ``End`` columns, which may fall on the
        day after ``Begin``.

    Returns
    -------
    `~astropy.time.Time` or `None`
        `None` is returned for missing (``////``) values.
    """
    raw = raw.strip()
    if not raw or "/" in raw:
        return None
    # A leading A (after), B (before) or U (uncertain) qualifier may precede
    # the HHMM value, e.g. "A0146" means the event began after 01:46 UT.
    digits = raw[-4:]
    hour, minute = int(digits[:2]), int(digits[2:])
    value = report_date + timedelta(hours=hour, minutes=minute)
    if after is not None and value < after:
        value += timedelta(days=1)
    return Time(value)


def _parse_event_lines(lines, report_date):
    """
    Parse the fixed-width event rows of a single day's report into a `~astropy.table.QTable`.
    """
    columns = {
        "event_id": [], "multiple_reports": [], "start_time": [], "max_time": [], "end_time": [],
        "observatory": [], "quality": [], "event_type": [], "location_or_frequency": [],
        "particulars": [], "region_number": [],
    }

    for line in lines:
        if _NO_EVENTS.search(line):
            break
        # Event rows are (at least) 80 characters wide. This also excludes
        # header/preamble lines that happen to start with digits, such as the
        # older file format's "HHMM UT DD Mon YYYY" creation timestamp line.
        if len(line) < 60 or not line[:4].strip().isdigit():
            continue

        raw = {name: line[start:end] for name, (start, end) in _COLUMNS.items()}

        region = raw["region_number"].strip()
        if region and not region.isdigit():
            # A very small number of historical rows (chiefly early CME
            # entries reported by the SOHO, "SOH", observatory) overflow the
            # standard 80 column width. Recover the region number from the
            # tail of the line instead.
            match = re.search(r"(\d+)\s*$", line)
            region = match.group(1) if match else ""

        start_time = _parse_time(raw["begin"], report_date)
        after = start_time.datetime if start_time else None
        columns["event_id"].append(int(raw["event_id"]))
        columns["multiple_reports"].append(raw["flag"].strip() == "+")
        columns["start_time"].append(start_time)
        columns["max_time"].append(_parse_time(raw["max"], report_date, after=after))
        columns["end_time"].append(_parse_time(raw["end"], report_date, after=after))
        columns["observatory"].append(raw["observatory"].strip())
        columns["quality"].append(raw["quality"].strip())
        columns["event_type"].append(raw["event_type"].strip())
        columns["location_or_frequency"].append(raw["location_or_frequency"].strip())
        columns["particulars"].append(" ".join(raw["particulars"].split()))
        columns["region_number"].append(int(region) if region else -1)

    table = QTable()
    # Explicit dtypes ensure that empty ("NO EVENT REPORTS.") tables still
    # have columns of the right type, so that they can be `vstack`-ed
    # together with non-empty tables in `parse_swpc_event_archive`.
    table["event_id"] = np.array(columns["event_id"], dtype=int)
    table["multiple_reports"] = np.array(columns["multiple_reports"], dtype=bool)
    table["start_time"] = Time(columns["start_time"]) if columns["start_time"] else Time([], format="isot")
    # `max_time` and `end_time` can be missing for a given row, so these are
    # kept as an object `Column` of `~astropy.time.Time` (or `None`)
    # instances, rather than letting `QTable` upcast them to a vectorized
    # `Time` mixin column (which it will only do when there happens to be no
    # `None` amongst the values, e.g. for an otherwise empty table).
    table["max_time"] = Column(columns["max_time"], dtype=object)
    table["end_time"] = Column(columns["end_time"], dtype=object)
    table["observatory"] = np.array(columns["observatory"], dtype=str)
    table["quality"] = np.array(columns["quality"], dtype=str)
    table["event_type"] = np.array(columns["event_type"], dtype=str)
    table["location_or_frequency"] = np.array(columns["location_or_frequency"], dtype=str)
    table["particulars"] = np.array(columns["particulars"], dtype=str)
    region_number = columns["region_number"]
    table["region_number"] = MaskedColumn(
        np.array(region_number, dtype=int), mask=[value == -1 for value in region_number]
    )
    table.meta["date"] = Time(report_date)
    return table


def _to_masked_time_column(column):
    """
    Convert an object `~astropy.table.Column` of `~astropy.time.Time` or
    `None` values (as produced by `_parse_event_lines`) into a masked
    ``Time`` column.

    This must only be applied as the very last step, after any
    `~astropy.table.vstack` calls: ``vstack`` silently discards the mask (and
    the underlying value) of an already-masked ``Time`` mixin column, rather
    than raising an error, so masking has to happen after all the tables
    that need combining have been.
    """
    if len(column) == 0:
        return Time([], format="isot")
    is_missing = np.array([value is None for value in column])
    filled = Time([value if value is not None else _MISSING_TIME_SENTINEL for value in column])
    filled[is_missing] = np.ma.masked
    return filled


def parse_swpc_event_file(filepath):
    """
    Parse a single NOAA SWPC "Edited Solar Events List" file.

    Parameters
    ----------
    filepath : `str` or `pathlib.Path`
        The path to a single day's ``events.txt`` file, e.g. as downloaded by
        `~sunkit_instruments.swpc.SWPCEventsClient`.

    Returns
    -------
    `~astropy.table.QTable`
        A table with one row per reported event.
        If the file reports no events for the day, an empty table is
        returned. In both cases, the UTC day the file describes is stored in
        ``.meta["date"]``.

    Notes
    -----
    * The ``Begin``, ``Max`` and ``End`` times are combined with the file's
      report date to construct full timestamps, following the day-rollover
      rules described in the file's README.
      Any leading ``A`` (after), ``B`` (before) or ``U`` (uncertain) qualifier
      on a time is dropped.
    * ``max_time`` and ``end_time`` are masked `~astropy.time.Time` columns,
      masked where the original report has a missing (``////``) value.
    * A very small number of historical reports (chiefly early CME entries
      reported by the SOHO, ``SOH``, observatory) exceed the standard
      80-column width. For these rows, ``region_number`` is recovered from
      the end of the line, but ``location_or_frequency`` and ``particulars``
      may be merged together.

    Examples
    --------
    >>> from sunkit_instruments.swpc import parse_swpc_event_file
    >>> parse_swpc_event_file(filepath)  # doctest: +SKIP
    """
    filepath = Path(filepath)
    lines = filepath.read_text(errors="replace").splitlines()
    report_date = _get_report_date(lines)
    table = _parse_event_lines(lines, report_date)
    # Safe to mask here, immediately, as this is the final table for a single
    # file - see the warning in `_to_masked_time_column` about combining
    # multiple tables with `~astropy.table.vstack` afterwards.
    table["max_time"] = _to_masked_time_column(table["max_time"])
    table["end_time"] = _to_masked_time_column(table["end_time"])
    table.meta["filename"] = filepath.name
    return table


def parse_swpc_event_archive(filepath):
    """
    Parse a yearly SWPC events archive from NOAA's warehouse.

    Parameters
    ----------
    filepath : `str` or `pathlib.Path`
        The path to a yearly ``{year}_events.tar.gz`` archive.

    Returns
    -------
    `~astropy.table.QTable`
        A table with one row per reported event across the whole year, sorted
        by ``start_time``.

    Examples
    --------
    >>> from sunkit_instruments.swpc import parse_swpc_event_archive
    >>> parse_swpc_event_archive(filepath)  # doctest: +SKIP
    """
    tables = []
    with tarfile.open(filepath, mode="r:gz") as tar:
        for member in sorted(tar.getmembers(), key=lambda m: m.name):
            if not member.isfile() or not re.match(r"^\d{8}events\.txt$", Path(member.name).name):
                continue
            lines = tar.extractfile(member).read().decode(errors="replace").splitlines()
            report_date = _get_report_date(lines)
            table = _parse_event_lines(lines, report_date)
            table.meta["filename"] = Path(member.name).name
            tables.append(table)

    if not tables:
        # Return an empty table with the expected schema.
        empty = QTable()
        empty["event_id"] = np.array([], dtype=int)
        empty["multiple_reports"] = np.array([], dtype=bool)
        empty["start_time"] = Time([], format="isot")
        empty["max_time"] = Time([], format="isot")
        empty["end_time"] = Time([], format="isot")
        empty["observatory"] = np.array([], dtype=str)
        empty["quality"] = np.array([], dtype=str)
        empty["event_type"] = np.array([], dtype=str)
        empty["location_or_frequency"] = np.array([], dtype=str)
        empty["particulars"] = np.array([], dtype=str)
        empty["region_number"] = MaskedColumn(np.array([], dtype=int), mask=[])
        return empty

    combined = vstack(tables, metadata_conflicts="silent")
    combined.sort("start_time")
    # Only convert `max_time`/`end_time` to masked `Time` columns now, after
    # all the vstacking above is done - see `_to_masked_time_column`.
    combined["max_time"] = _to_masked_time_column(combined["max_time"])
    combined["end_time"] = _to_masked_time_column(combined["end_time"])
    return combined


def parse_swpc_events(filepath, time_range=None):
    """
    Parse a NOAA SWPC events file, whichever kind it is.

    This dispatches automatically to `parse_swpc_event_file` or
    `parse_swpc_event_archive`, depending on whether ``filepath`` is a single
    day's ``events.txt`` file or a yearly ``{year}_events.tar.gz`` archive.
    This is the recommended way to parse files downloaded with
    `~sunkit_instruments.swpc.SWPCEventsClient`: depending on how old the
    requested data is, a search may be satisfied by either kind of file, and
    this avoids having to inspect the downloaded file yourself to know which
    parser to call.

    Parameters
    ----------
    filepath : `str` or `pathlib.Path`
        The path to a SWPC events file, either a single day's ``events.txt``
        file or a yearly ``{year}_events.tar.gz`` archive.
    time_range : `~sunpy.time.TimeRange`, optional
        If given, only events with a ``start_time`` inside this range are
        returned. This is particularly useful when ``filepath`` is a yearly
        archive: a single search result for an archived year always
        downloads, and therefore parses, the whole year, which may be far
        more than was originally asked for.

    Returns
    -------
    `~astropy.table.QTable`

    Examples
    --------
    >>> from sunkit_instruments.swpc import parse_swpc_events
    >>> parse_swpc_events(filepath)  # doctest: +SKIP
    """
    filepath = Path(filepath)
    table = parse_swpc_event_archive(filepath) if tarfile.is_tarfile(filepath) else parse_swpc_event_file(filepath)
    if time_range is not None:
        table = table[(table["start_time"] >= time_range.start) & (table["start_time"] <= time_range.end)]
    return table
