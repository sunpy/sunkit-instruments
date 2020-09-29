"""
Sunkit-instruments
==================

A SunPy-affiliated package for solar instrument-specific tools.

* Homepage: https://sunpy.org
* Documentation: https://sunkit-instruments.readthedocs.io/en/latest/
"""
import sys

# Enforce Python version check during package import.
__minimum_python_version__ = "3.7"

try:
    from .version import __version__
except ImportError:
    print("version.py not found, please reinstall sunkit_instruments.")
    __version__ = "unknown"


class UnsupportedPythonError(Exception):
    """
    Running on an unsupported version of Python.
    """


if sys.version_info < tuple(int(val) for val in __minimum_python_version__.split('.')):
    # This has to be .format to keep backwards compatibly.
    raise UnsupportedPythonError(
        "sunkit_instruments does not support Python < {}".format(__minimum_python_version__))

__all__ = []
