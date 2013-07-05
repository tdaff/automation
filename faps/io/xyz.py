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

class XYZReader(StructureReader):
    """
    XYZ file reader with some symmetry and custom properties.

    """

    filetypes = ['*.cif', 'cif']

    @staticmethod
    def read_file(filename=None, string=None, template=None):
        """Read a structure from an xyz file."""
        info("Reading positions from xyz file: %s" % filename)
        filetemp = open(filename)
        xyz_file = filetemp.readlines()
        filetemp.close()
        # Setting the cell
        if len(cell) == 6:
            self.cell.params = cell
        elif len(cell) == 9:
            self.cell.cell = array(cell).reshape((3, 3))
        elif cell is not None:
            error("Invalid cell specification %s" % str(cell))
        # Build a local list before setting attribute
        newatoms = []
        natoms = int(xyz_file[0])
        self.properties['header'] = xyz_file[1].strip()
        for line in xyz_file[2:2+natoms]:
            newatom = Atom()
            newatom.from_xyz(line)
            newatoms.append(newatom)
        if update:
            if natoms != self.natoms:
                critical("Incorrect number of atoms to update")
                terminate(96)
            for atom, newatom in zip(self.atoms, newatoms):
                if atom.type != newatom.type:
                    error("Atom order may have changed")
                atom.pos = newatom.pos
        else:
            self.atoms = newatoms
            self.order_by_types()

        return new_structure