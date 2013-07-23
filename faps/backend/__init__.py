"""
faps.backend module contains interfaces for loading and saving
simulation states with files or databases.

"""

from abc import ABCMeta


class Backend(object):
    """
    Simple abstract base class for backends.

    """

    __metaclass__ = ABCMeta

    filetypes = {}

    @staticmethod
    @abstractmethod
    def load(jobname):
        return None

    @staticmethod
    @abstractmethod
    def save(jobname):
        return None
