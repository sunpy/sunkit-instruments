__all__ = ["AbstractChannel"]

import abc


class AbstractChannel(abc.ABC):

    @abc.abstractmethod
    def temperature_response(self):
        ...

    @abc.abstractmethod
    def effective_area(self):
        ...
