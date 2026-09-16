from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ["SWPCEventsClient"]


class SWPCEventsClient(GenericClient):
    """
    Provides access to the NOAA SWPC "Edited Solar Events List".

    These are daily lists of solar flares, radio bursts and other transient
    solar-terrestrial events, compiled from preliminary reports received at
    the Space Weather Prediction Center (SWPC) and manually reviewed and
    edited by the duty forecaster.
    See the `README <ftp://ftp.swpc.noaa.gov/pub/indices/events/README>`__
    for a full description of the data and its history.

    Data is available from 1966 onwards, but machine readable data is only
    available from 1996 onwards, gathered from three locations on
    `SWPC's FTP server <ftp://ftp.swpc.noaa.gov/pub/>`__:

    * complete previous years, bundled into a single
      ``pub/warehouse/{year}/{year}_events.tar.gz`` archive per year;
    * the current (not yet complete) year, from
      ``pub/warehouse/{year}/{year}_events/``;
    * the most recent files, from ``pub/indices/events/``. The `README
      <ftp://ftp.swpc.noaa.gov/pub/indices/events/README>`__ describes this as
      a rolling ~60 day window, though in practice a much longer history is
      typically retained.

    The yearly archive is preferred whenever a requested year already has
    one - a single download beats scraping and fetching potentially
    hundreds of individual daily files for the same year - falling back to
    daily files only for a year that isn't archived yet (in practice, just
    the current year). A single query can therefore return a mix of yearly
    archives and single-day files - use
    `sunkit_instruments.swpc.parse_swpc_events` to parse either kind without
    having to check which one you got. Note that a result for an archived
    year always downloads (and therefore parses) the whole year: pass
    ``time_range`` to `~sunkit_instruments.swpc.parse_swpc_events` if you
    only want the events from your original, narrower search.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> import sunkit_instruments.swpc  # doctest: +SKIP
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.swpc_events)  # doctest: +SKIP
    >>> results  # doctest: +SKIP
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the SWPCEventsClient:
    Source: ftp://ftp.swpc.noaa.gov/pub/indices/events/
    <BLANKLINE>
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ----------- ... ------ --------
    2016-01-01 00:00:00.000 2016-12-31 23:59:59.999 SWPC-EVENTS ...   SWPC     NOAA
    <BLANKLINE>
    <BLANKLINE>
    >>> results[0]["url"]  # doctest: +SKIP
    'ftp://ftp.swpc.noaa.gov/pub/warehouse/2016/2016_events.tar.gz'
    >>> from sunkit_instruments.swpc import parse_swpc_events
    >>> files = Fido.fetch(results)  # doctest: +SKIP
    >>> events = [parse_swpc_events(file) for file in files]  # doctest: +SKIP
    """
    # The most recent files, refreshed multiple times a day.
    pattern_recent = "ftp://ftp.swpc.noaa.gov/pub/indices/events/{{year:4d}}{{month:2d}}{{day:2d}}events.txt"
    # The current (not yet archived) year, one file per day.
    pattern_current_year = ("ftp://ftp.swpc.noaa.gov/pub/warehouse/{{year:4d}}/{{year:4d}}_events/"
                             "{{year:4d}}{{month:2d}}{{day:2d}}events.txt")
    # Complete previous years, bundled as a single yearly archive.
    pattern_archive = "ftp://ftp.swpc.noaa.gov/pub/warehouse/{{year:4d}}/{{year:4d}}_events.tar.gz"

    @property
    def info_url(self):
        return "ftp://ftp.swpc.noaa.gov/pub/indices/events/"

    def search(self, *args, **kwargs):
        """
        Query this client for a list of results.

        Parameters
        ----------
        \\*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        \\*\\*kwargs: `dict`
             Any extra keywords to refine the search.

        Returns
        -------
        `~sunpy.net.dataretriever.QueryResponse`
            The query result, combining daily files and yearly archives as
            described above.
        """
        matchdict = self._get_match_dict(*args, **kwargs)
        tr = TimeRange(matchdict["Start Time"], matchdict["End Time"])
        metalist = []

        # A completed year's archive is a single file covering the whole
        # year, whereas the daily-file patterns below can mean hundreds of
        # individual FTP downloads for the same year - since `pub/indices/events/`
        # in practice retains far more than the ~60 days the README
        # documents (see the class docstring). So the archive is preferred
        # whenever a requested year already has one.
        seen_years = set()
        scraper = Scraper(format=self.pattern_archive)
        for exdict in scraper._extract_files_meta(tr):
            seen_years.add(exdict["year"])
            metalist.append(self.post_search_hook(exdict, matchdict))

        # Any day in a year with no archive yet - in practice, just the
        # current, not-yet-complete year - falls back to the daily files.
        # The same day can be listed both in the rolling recent window and
        # in the current year's warehouse directory, so results found via
        # `pattern_recent` take priority and duplicates are dropped.
        seen_days = set()
        for pattern in (self.pattern_recent, self.pattern_current_year):
            scraper = Scraper(format=pattern)
            for exdict in scraper._extract_files_meta(tr):
                if exdict["year"] in seen_years:
                    continue
                key = (exdict["year"], exdict["month"], exdict["day"])
                if key in seen_days:
                    continue
                seen_days.add(key)
                metalist.append(self.post_search_hook(exdict, matchdict))

        metalist = sorted(metalist, key=lambda row: row["url"])
        return QueryResponse(metalist, client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {
            attrs.Instrument: [
                ("SWPC-EVENTS", "NOAA SWPC Edited Solar Events List.")],
            attrs.Physobs: [
                ("event_list", "A catalog of solar flares, radio bursts, and other transient solar-terrestrial events.")],
            attrs.Source: [("SWPC", "The Space Weather Prediction Center.")],
            attrs.Provider: [("NOAA", "The National Oceanic and Atmospheric Administration.")],
        }
        return adict
