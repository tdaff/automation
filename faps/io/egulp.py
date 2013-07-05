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
from . import StructureReader, StructureWriter


class EgulpReader(StructureReader):
    def charges_from_egulp(self, filename):
        """Parse QEq charges from EGULP output."""
        info("Getting charges from file: %s" % filename)
        filetemp = open(filename)
        gout = filetemp.readlines()
        filetemp.close()
        # charges.dat file only has list of charges in it
        for atom, chg_line in zip(self.atoms, gout):
            atom.charge = float(chg_line.split()[2])
            if atom.charge != atom.charge:
                # These are 'nan'; should not continue or fastmc will mess up
                error("Egulp charges did not converge, check structure")
                terminate(107)
            elif abs(atom.charge) == INFINITY:
                error("Egulp gave infinite charges, check structure")
                terminate(108)
            elif abs(atom.charge) > 10:
                warning("Very high charge from egulp: %s %f"
                        (atom.site, atom.charge))




class EgulpWriter(StructureWriter):
    """
    CIF file reader with some symmetry and custom properties.

    """

    filetypes = ['*.cif', 'cif']

    @staticmethod
    def to_egulp(self, typed_atoms=False):
        """Generate input files for Eugene's QEq code."""
        # bind cell locally for speed and convenience
        cell = self.cell.cell
        inv_cell = self.cell.inverse
        # format exactly as Eugene original script generates it
        geometry_file = ['%s\n' % self.name]
        geometry_file.extend(self.cell.to_vector_strings(fmt=' %15.12f'))
        geometry_file.append('%i\n' % self.natoms)

        atomic_numbers = self.atomic_numbers

        if typed_atoms:
            # Include custom typing, new types must be found by hand.
            # 800 N=O -- removed
            # 801 S=O
            # 802 S-O-H
            for atom_idx, atom in enumerate(self.atoms):
                if atom.uff_type == 'S_3+6':
                    for bond in self.bonds:
                        if atom_idx in bond:
                            other_idx = other_bond_index(bond, atom_idx)
                            if self.atoms[other_idx].uff_type == 'O_2':
                                atomic_numbers[other_idx] = 801
                            elif self.atoms[other_idx].uff_type == 'O_3':
                                atomic_numbers[other_idx] = 802
                                for another_bond in self.bonds:
                                    if other_idx in another_bond:
                                        another_idx = other_bond_index(another_bond, other_idx)
                                        if self.atoms[another_idx].uff_type == 'H_':
                                            atomic_numbers[another_idx] = 1001
# TODO(tdaff): delete code in future version; removed NO2 typing
#                elif atom.uff_type == 'N_R':
#                    this_bonds = []
#                    for bond in self.bonds:
#                        if atom_idx in bond:
#                            other_idx = other_bond_index(bond, atom_idx)
#                            if self.atoms[other_idx].uff_type == 'O_R':
#                                atomic_numbers[other_idx] = 800

        # atomic numbers should have been modified with exotic types by now
        geometry_file.extend([
            ('%6d ' % atomic_number) +
            ('%12.7f %12.7f  %12.7f' % tuple(atom.ipos(cell, inv_cell))) +
            ('%12.7f\n' % atom.charge)
            for atom, atomic_number in zip(self.atoms, atomic_numbers)
        ])

        return geometry_file
