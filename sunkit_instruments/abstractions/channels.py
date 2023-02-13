__all__ = ["AbstractChannel"]

import abc

import astropy.units as u


class AbstractChannel(abc.ABC):
    @abc.abstractmethod
    def effective_area(self) -> u.cm**2:
        ...

    @abc.abstractmethod
    def gain(self) -> u.DN / u.photon:
        ...

    @abc.abstractmethod
    def temperature_response(self) -> u.cm**5 * u.ct / (u.pix * u.s):
        ...
