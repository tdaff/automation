"""

Abstract classes that can be plugged in and register filetypes.

"""

from abc import abstractmethod
from abc import ABCMeta


class StructureReader(object):
    """
    Abstract base class for file readers.
    Plugins must be a subclass of StructureReader and
    register their filetypes.

    """

    __metaclass__ = ABCMeta

    filetypes = {}

    def __init__(self, ):
        self.register_filetypes()

    @classmethod
    def register_filetypes(cls):
        for subclass in cls.__subclasses__():
            for filetype in subclass.filetypes:
                cls.filetypes[filetype] = subclass

    @staticmethod
    @abstractmethod
    def read_file(filename):
        return None

    @staticmethod
    @abstractmethod
    def update_from_file(structure, filename):
        return None


class StructureWriter(object):
    """
    Abstract base class for structure file writers.
    Plugins must be a subclass of StructureWriter and
    register their filetypes.

    """

    __metaclass__ = ABCMeta

    filetypes = {}

    def __init__(self, ):
        self.register_filetypes()

    @classmethod
    def register_filetypes(cls):
        for subclass in cls.__subclasses__():
            for filetype in subclass.filetypes:
                cls.filetypes[filetype] = subclass

    @staticmethod
    @abstractmethod
    def write_file(structure, filename):
        return None
