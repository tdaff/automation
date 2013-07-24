"""
Register arbritary simulation methods and processors by subclassing
the abstract method.

"""

from abc import abstractmethod
from abc import ABCMeta


class Calculator(object):
    """
    Abstract base class for calculators.
    Calculators must be a subclass of Calculator and
    register the attributes that they 'provides'.

    """

    __metaclass__ = ABCMeta

    provides = {}
    requires = {}

    def __init__(self):
        self.register_calculators()

    @classmethod
    def register_calculators(cls):
        for subclass in cls.__subclasses__():
            for provide in subclass.provides:
                cls.provides[filetype] = subclass
            for requires in subclass.requires:
                cls.requires[filetype] = subclass

    @staticmethod
    @abstractmethod
    def run_calculation(self):
        return None
