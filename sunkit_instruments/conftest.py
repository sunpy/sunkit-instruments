import pytest


@pytest.fixture()
def sunpy_cache(mocker, tmp_path):
    """
    Provide a way to add local files to the cache.

    This can be useful when mocking remote requests.
    """
    from types import MethodType

    from sunpy.data.data_manager.cache import Cache
    from sunpy.data.data_manager.downloader import ParfiveDownloader
    from sunpy.data.data_manager.storage import InMemStorage

    cache = Cache(ParfiveDownloader(), InMemStorage(), tmp_path, None)

    def add(self, url, path):
        self._storage.store(
            {
                "url": url,
                "file_path": path,
                "file_hash": "none",  # hash doesn't matter
            }
        )

    cache.add = MethodType(add, cache)

    def func(mocked):
        mocker.patch(mocked, cache)
        return cache

    yield func
