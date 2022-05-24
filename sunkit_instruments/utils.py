"""
This module provides general utility functions.
"""
import shutil
import tarfile
from pathlib import Path

import requests


def _download_data(urls: list):
    """
    This only should be called in the example gallery.

    THIS IS NOT A USER FACING FUNCTION TO DOWNLOAD DATA.
    """
    for url in urls:
        filename = url.split("/")[-1]
        if not Path(filename).exists():
            with requests.get(url, stream=True) as r:
                with open(filename, "wb") as f:
                    shutil.copyfileobj(r.raw, f)
        if ".tar.gz" in filename:
            with tarfile.open(filename, "r") as tar:
                tar.extractall(Path(filename).parent)
