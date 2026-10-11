from astropy.time import Time

from sunpy.time import TimeRange

from sunkit_instruments.data.test import get_test_filepath
from sunkit_instruments.swpc import parse_swpc_event_archive, parse_swpc_event_file, parse_swpc_events

NEW_FORMAT_FILE = get_test_filepath("swpc_20260710events.txt")
OLD_FORMAT_FILE = get_test_filepath("swpc_19970121events.txt")
NO_EVENTS_FILE = get_test_filepath("swpc_19970115events.txt")
ARCHIVE_FILE = get_test_filepath("swpc_1997_events_sample.tar.gz")


def test_parse_new_format_file():
    table = parse_swpc_event_file(NEW_FORMAT_FILE)
    assert table.meta["date"] == Time("2026-07-10")
    assert len(table) == 18
    assert list(table["event_id"][:4]) == [7120, 7120, 7120, 7120]
    assert table["multiple_reports"][2]
    assert not table["multiple_reports"][0]
    assert table["observatory"][0] == "PAL"
    assert table["event_type"][0] == "RSP"
    assert table["location_or_frequency"][0] == "025-075"
    assert table["particulars"][0] == "VI/2"
    # No region number was reported for the first event.
    assert table["region_number"].mask[0]


def test_parse_new_format_multi_value_particulars():
    table = parse_swpc_event_file(NEW_FORMAT_FILE)
    row = table[table["event_id"] == 7170][0]
    assert row["event_type"] == "XRA"
    assert row["particulars"] == "C6.0 1.1E-02"
    assert row["region_number"] == 4485


def test_parse_day_rollover():
    table = parse_swpc_event_file(NEW_FORMAT_FILE)
    row = table[table["event_id"] == 7230][0]
    # Begin is 11:10 and End is 00:10, so End must roll over to the next day.
    assert row["start_time"] == Time("2026-07-10T11:10:00")
    assert row["end_time"] == Time("2026-07-11T00:10:00")


def test_parse_missing_max_time():
    table = parse_swpc_event_file(NEW_FORMAT_FILE)
    assert isinstance(table["max_time"], Time)
    assert table["max_time"].mask[0]


def test_parse_old_format_file():
    table = parse_swpc_event_file(OLD_FORMAT_FILE)
    assert table.meta["date"] == Time("1997-01-21")
    assert len(table) == 6
    assert table["event_id"][0] == 5060
    assert table["start_time"][0] == Time("1997-01-21T00:21:00")
    assert table["end_time"][0] == Time("1997-01-21T00:21:00")
    assert table["max_time"].mask.all()
    assert table["region_number"].mask.all()


def test_parse_no_events_file():
    table = parse_swpc_event_file(NO_EVENTS_FILE)
    assert table.meta["date"] == Time("1997-01-15")
    assert len(table) == 0


def test_parse_swpc_event_archive():
    table = parse_swpc_event_archive(ARCHIVE_FILE)
    # 6 events on the 21st, and 0 on the 15th.
    assert len(table) == 6
    assert set(table["event_id"]) == {5060, 5070, 5080, 5090, 5100, 5120}


def test_parse_swpc_event_archive_masks_missing_times():
    # `parse_swpc_event_archive` combines multiple files with
    # `~astropy.table.vstack`, so its masking of `max_time`/`end_time` has to
    # happen after that `vstack` call, since `vstack` does not preserve the
    # mask of an already-masked `Time` column.
    table = parse_swpc_event_archive(ARCHIVE_FILE)
    assert isinstance(table["max_time"], Time)
    assert isinstance(table["end_time"], Time)
    # All 6 events on the 21st are missing a Max time, but have a real End time.
    assert table["max_time"].mask.all()
    assert not table["end_time"].mask.any()
    assert table["end_time"][0] == Time("1997-01-21T00:21:00")


def test_parse_swpc_events_dispatches_on_file_type():
    file_table = parse_swpc_events(OLD_FORMAT_FILE)
    assert len(file_table) == 6

    archive_table = parse_swpc_events(ARCHIVE_FILE)
    assert len(archive_table) == 6


def test_parse_swpc_events_time_range():
    # Restrict the archive (both days) down to just the 6 events on the 21st.
    time_range = TimeRange("1997-01-20", "1997-01-22")
    table = parse_swpc_events(ARCHIVE_FILE, time_range=time_range)
    assert len(table) == 6
    assert set(table["event_id"]) == {5060, 5070, 5080, 5090, 5100, 5120}

    # A range that only covers the (event-free) 15th should return nothing.
    time_range = TimeRange("1997-01-15", "1997-01-16")
    table = parse_swpc_events(ARCHIVE_FILE, time_range=time_range)
    assert len(table) == 0
