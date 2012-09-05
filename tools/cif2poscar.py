#!/usr/bin/env python

"""
cif2poscar.py

Convert a cif file to a vasp5 poscar

"""

import os
import sys
import time

# faps is in the parent directory of this script
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from faps import Structure, unique


def to_vasp(structure):
    """Return a vasp5 poscar as a list of lines."""
    types = [atom.type for atom in structure.atoms]
    ordered_types = unique(types)
    poscar = ["%s\n" % structure.name[:80],
              " 1.0\n"]
    # Vasp does 16 dp but we get rounding errors (eg cubic) use 14
    poscar.extend(structure.cell.to_vector_strings(fmt="%23.14f"))
    poscar.append("".join("%5s" % x for x in ordered_types) + "\n")
    poscar.append("".join("%6i" % types.count(x)
                          for x in ordered_types) + "\n")
    poscar.extend(["Cartesian\n"])

    # assume atoms are ordered
    # Atom lines differ from vasp to ensure spaces between numbers
    # with <-10 values
    for atom in structure.atoms:
        poscar.append("%20.15f %19.15f %19.15f\n" % tuple(atom.pos))

    return poscar


def main():
    """
    Run the conversion.

    """
    filename = sys.argv[1]
    cif_file = Structure("cif file")
    cif_file.from_cif(filename)
    output_file = open(sys.argv[2], 'w')
    output_file.writelines(to_vasp(cif_file))


if __name__ == '__main__':

    main()
