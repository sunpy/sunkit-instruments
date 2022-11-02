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
                
                import os
                
                def is_within_directory(directory, target):
                    
                    abs_directory = os.path.abspath(directory)
                    abs_target = os.path.abspath(target)
                
                    prefix = os.path.commonprefix([abs_directory, abs_target])
                    
                    return prefix == abs_directory
                
                def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                
                    for member in tar.getmembers():
                        member_path = os.path.join(path, member.name)
                        if not is_within_directory(path, member_path):
                            raise Exception("Attempted Path Traversal in Tar File")
                
                    tar.extractall(path, members, numeric_owner=numeric_owner) 
                    
                
                safe_extract(tar, Path(filename).parent)
