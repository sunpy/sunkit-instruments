"""This module defines abstractions related to instrument channels."""

__all__ = ["AbstractChannel"]

import abc

import astropy.units as u


class AbstractChannel(abc.ABC):
    """An abstract base class for defining instrument channels.

    .. caution::

       This abstract class is still under development and may change
       in the near future.
    """

    @abc.abstractmethod
    def effective_area(self) -> u.cm**2:
        ...

    @abc.abstractmethod
    def gain(self) -> u.DN / u.photon:
        ...

    @abc.abstractmethod
    def temperature_response(self) -> u.cm**5 * u.ct / (u.pix * u.s):
        ...
