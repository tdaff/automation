#!/usr/bin/env python

"""
function_switch.py

Alter a known structure with new functional groups ready for fapping.

"""

import os
import sys
import time

# faps is in the parent directory of this script
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from faps import Structure


def to_cif(structure):
    """Return a CIF file without bonding or atom types."""

    atoms = structure.atoms
    cell = structure.cell
    name = "converted structure"
    inv_cell = cell.inverse

    type_count = {}

    atom_part = []
    for atom in atoms:
        if atom is None:
            # blanks are left in here
            continue
        if hasattr(atom, 'uff_type') and atom.uff_type is not None:
            uff_type = atom.uff_type
        else:
            uff_type = '?'
        if atom.element in type_count:
            type_count[atom.element] += 1
        else:
            type_count[atom.element] = 1
        atom.site = "%s%i" % (atom.element, type_count[atom.element])
        atom_part.append("%-5s %-5s %-5s " % (atom.site, atom.element, uff_type))
        atom_part.append("%f %f %f\n" % tuple(atom.ifpos(inv_cell)))

    cif_file = [
        "data_%s\n" % name.replace(' ', '_'),
        "%-33s %s\n" % ("_audit_creation_date", time.strftime('%Y-%m-%dT%H:%M:%S%z')),
        "%-33s %s\n" % ("_audit_creation_method", "poscar2cif"),
        "%-33s %s\n" % ("_symmetry_space_group_name_H-M", "P1"),
        "%-33s %s\n" % ("_symmetry_Int_Tables_number", "1"),
        "%-33s %s\n" % ("_space_group_crystal_system", cell.crystal_system),
        "%-33s %f\n" % ("_cell_length_a", cell.a),
        "%-33s %f\n" % ("_cell_length_b", cell.b),
        "%-33s %f\n" % ("_cell_length_c", cell.c),
        "%-33s %f\n" % ("_cell_angle_alpha", cell.alpha),
        "%-33s %f\n" % ("_cell_angle_beta", cell.beta),
        "%-33s %f\n" % ("_cell_angle_gamma", cell.gamma),
        # start of atom loops
        "\nloop_\n",
        "_atom_site_label\n",
        "_atom_site_type_symbol\n",
        "_atom_type_description\n",
        "_atom_site_fract_x\n",
        "_atom_site_fract_y\n",
        "_atom_site_fract_z\n"] + atom_part

    return cif_file


def main():
    """
    Run the conversion.

    """
    filename = sys.argv[1]
    poscar = Structure("vasp structure")
    poscar.from_vasp(filename)
    output_file = open(sys.argv[2], 'w')
    output_file.writelines(to_cif(poscar))


if __name__ == '__main__':

    main()
