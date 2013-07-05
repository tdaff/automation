"""

Reader for VASP coordinate files (Crystallographic information files).

Parses symmetry and some custom properties.

"""
import shlex
from logging import debug, info, warning, error

from faps.core.structure import Structure
from faps.core.symmetry import Symmetry
from faps.core.cell import Cell
from faps.util.text import strip_blanks, ufloat
from . import StructureReader

class PoscarReader(StructureReader):
    """
    CIF file reader with some symmetry and custom properties.

    """

    filetypes = ['*.cif', 'cif']

    @staticmethod
    def update_from_file(structure, filename, symmetry=False):
        """Parse charges and update structure."""
        info("Getting charges from file: %s" % filename)
        charges = []
        filetemp = open(filename)
        for line in filetemp:
            if line.startswith(" Charge"):
                line = line.split()
                # index, type, charge
                charges.append((int(line[1]), int(line[4]), float(line[6])))
            if "Error" in line:
                if float(line.split()[-1]) > 0.6:
                    warning("Error in repeat charges is very high - check cube!")
        filetemp.close()
        if symmetry:
            tree = self.symmetry_tree
            if len(charges) != len(tree):
                error("Incorrect number of charge sets; check REPEAT output")
                terminate(97)
            for symm, charge in zip(sorted(tree.items()), charges):
                for at_idx in symm[1]:
                    self.atoms[at_idx].charge = charge[2]
        else:
            if len(charges) != len(self.atoms):
                error("Incorrect number of charges; check REPEAT output")
                terminate(90)
            for atom, charge in zip(self.atoms, charges):
                atom.charge = charge[2]
