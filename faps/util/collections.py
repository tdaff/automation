"""
Collection of useful utility objects

Classes with general use cases are included here, similar to the
sdlib collections module.

"""

class AttrDict(dict):
    """Dictionary object that allows attribute based lookup."""
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
   
