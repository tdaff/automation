#!/usr/bin/env python

"""
Guest object for setting up Monte Carlo simulations.

"""

try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import os
import sys
from logging import warning, debug, error, info, critical
from os import path

import numpy as np

from elements import WEIGHT, ATOMIC_NUMBER, UFF, VASP_PSEUDO_PREF
from elements import CCDC_BOND_ORDERS, GULP_BOND_ORDERS, METALS
from elements import COVALENT_RADII, UFF_FULL, QEQ_PARAMS

from faps.util.text import strip_blanks


class Guest(object):
    """Guest molecule and properties."""
    def __init__(self, ident=None, guest_path=None):
        """Populate an empty guest then load from library if required."""
        self.ident = ''
        self.name = "Unknown guest"
        self.potentials = {}
        self.probability = []
        self.atoms = []
        self.source = "Unknown source"
        self.uptake = {}
        self.hoa = {}
        self.c_v = {}
        # only load if asked, set the ident in the loader
        if ident:
            self.load_guest(ident, guest_path=None)

    def load_guest(self, ident, guest_path=None):
        """Look in guests.lib in submit directory and default."""
        # Ident set here to keep consistent
        self.ident = ident
        # Need the different directories
        if guest_path is None:
            guest_path = os.getcwd()
        if __name__ != '__main__':
            script_dir = path.dirname(__file__)
        else:
            script_dir = path.abspath(sys.path[0])
        dot_faps_dir = path.join(path.expanduser('~'), '.faps')
        # A parser for each location
        job_guests = configparser.SafeConfigParser()
        dot_faps_guests = configparser.SafeConfigParser()
        lib_guests = configparser.SafeConfigParser()
        # Try and find guest in guests.lib
        job_guests.read(path.join(guest_path, 'guests.lib'))
        dot_faps_guests.read(path.join(dot_faps_dir, 'guests.lib'))
        lib_guests.read(path.join(script_dir, 'guests.lib'))
        # Job dir and user defined have priority
        if job_guests.has_section(ident):
            debug("%s found in job dir" % ident)
            self._parse_guest(job_guests.items(ident))
        elif dot_faps_guests.has_section(ident):
            debug("%s found in dot_faps" % ident)
            self._parse_guest(dot_faps_guests.items(ident))
        elif lib_guests.has_section(ident):
            debug("%s found in library" % ident)
            self._parse_guest(lib_guests.items(ident))
        else:
            error("Guest not found: %s" % ident)

    def _parse_guest(self, raw_text):
        """Set attributes according to the raw input."""
        for key, val in raw_text:
            if key == 'atoms':
                # Build a local list to replace what is there
                new_atoms = []
                # Only use non blank lines
                atoms = strip_blanks(val.splitlines())
                for atom in atoms:
                    atom = atom.split()
                    new_atoms.append(Atom(
                        at_type=atom[0],
                        mass=float(atom[1]),
                        charge=float(atom[2]),
                        pos=tuple(float(x) for x in atom[3:6])))
                self.atoms = new_atoms
            elif key == 'potentials':
                potens = strip_blanks(val.splitlines())
                for poten in potens:
                    poten = poten.split()
                    self.potentials[poten[0]] = tuple(float(x)
                                                      for x in poten[1:])
            elif key == 'probability':
                prob = strip_blanks(val.splitlines())
                self.probability = [tuple(int(y) for y in x.split())
                                    for x in prob]
            elif key == 'molar volume':
                self._molar_volume = float(val)
            elif key == 'probe radius':
                self.probe_radius = float(val)
            else:
                # Arbitrary attributes can be set
                # Might also enable breakage
                setattr(self, key, val)

    # Simple attribute-like calculations
    @property
    def types(self):
        """Ordered list of atom types."""
        return [atom.type for atom in self.atoms]

    @property
    def weight(self):
        """Unit cell weight."""
        return sum([atom.mass for atom in self.atoms])

    @property
    def molar_volume(self):
        """Molar volume at STP in dm3/mol"""
        if hasattr(self, '_molar_volume'):
            return self._molar_volume
        else:
            return 22.414

    @property
    def to_dict(self):
        drepr = {"@module": self.__class__.__module__,
                 "@class": self.__class__.__name__,
                 "cell": self.cell.tolist()}
        return drepr

    @staticmethod
    def dict(drepr):
        return Cell(cell=drepr['cell'])

