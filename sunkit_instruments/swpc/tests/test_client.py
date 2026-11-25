from datetime import date, timedelta

import pytest

from sunpy.net import Fido
from sunpy.net import attrs as a

from sunkit_instruments.swpc import SWPCEventsClient


def test_can_handle_query():
    assert SWPCEventsClient._can_handle_query(
        a.Time("2016/1/1", "2016/1/2"), a.Instrument.swpc_events
    )
    assert not SWPCEventsClient._can_handle_query(a.Time("2016/1/1", "2016/1/2"))
    assert not SWPCEventsClient._can_handle_query(a.Instrument.swpc_events)


@pytest.mark.remote_data
def test_fido_search_daily_files():
    # Use recent dates (rather than a hardcoded year) so this always lands in
    # the current, not-yet-archived year and exercises the daily-file path.
    # A few days in the past, rather than today/yesterday, to stay clear of
    # any short delay in that day's file being published.
    end = date.today() - timedelta(days=2)
    start = end - timedelta(days=1)
    results = Fido.search(
        a.Time(str(start), str(end)), a.Instrument.swpc_events
    )
    client_results = results[0]
    assert len(client_results) == 2
    assert set(client_results["Instrument"]) == {"SWPC-EVENTS"}
    assert set(client_results["Provider"]) == {"NOAA"}


@pytest.mark.remote_data
def test_fido_search_archived_year():
    # 1997 is long since archived into a single yearly tar.gz.
    results = Fido.search(
        a.Time("1997/1/1", "1997/1/2"), a.Instrument.swpc_events
    )
    client_results = results[0]
    assert len(client_results) == 1
    assert client_results[0]["url"].endswith("1997_events.tar.gz")


@pytest.mark.remote_data
def test_fido_search_prefers_archive_over_daily_files():
    # 2016 is complete and archived, even though its individual daily files
    # also still happen to exist in the rolling "recent" directory (which in
    # practice retains far more than the ~60 days the README describes) -
    # the single yearly archive should still be preferred, rather than the
    # much less efficient one-row-per-day alternative.
    results = Fido.search(
        a.Time("2016/1/1", "2016/1/2"), a.Instrument.swpc_events
    )
    client_results = results[0]
    assert len(client_results) == 1
    assert client_results[0]["url"].endswith("2016_events.tar.gz")


@pytest.mark.remote_data
def test_fido_fetch_recent(tmp_path):
    end = date.today() - timedelta(days=2)
    results = Fido.search(
        a.Time(str(end), str(end)), a.Instrument.swpc_events
    )
    files = Fido.fetch(results[0], path=tmp_path / "{file}")
    assert len(files) == 1
