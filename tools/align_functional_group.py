"""
create_functional_group.py

Generate a functional group block for a functional_groups.lib file. Takes a
smiles string and uses openbabel to deduce the 3D structure and aligns it
based on a benzene ring.

"""

import argparse

import openbabel as ob
import numpy as np
import pybel
from numpy import asarray, cross, dot, array
from numpy.linalg import norm


ATOMIC_NUMBER = [
    "ZERO", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uuo"]


def matrix_rotate(source, target):
    """Create a rotation matrix that will rotate source on to target."""
    # Normalise so there is no scaling in the array
    source = asarray(source)/norm(source)
    target = asarray(target)/norm(target)
    v = cross(source, target)
    vlen = dot(v, v)
    if vlen == 0.0:
        # already aligned, no rotation needed
        return identity(3)
    c = dot(source, target)
    h = (1 - c)/(np.dot(v, v))
    return array([[c + h*v[0]*v[0], h*v[0]*v[1] - v[2], h*v[0]*v[2] + v[1]],
                  [h*v[0]*v[1] + v[2], c + h*v[1]*v[1], h*v[1]*v[2] - v[0]],
                  [h*v[0]*v[2] - v[1], h*v[1]*v[2] + v[0], c + h*v[2]*v[2]]])


def realign(coordinates, origin_index, up_index, plane_index):
    # move everything
    array_coords = asarray(coordinates)
    array_coords -= array_coords[origin_index]
    # align with axis
    current_axis = array_coords[up_index] - array_coords[origin_index]
    align_matrix = matrix_rotate(-current_axis, [0, 1, 0])
    array_coords = asarray([dot(align_matrix, x) for x in array_coords])
    # align with plane
    in_plane_vector = array_coords[plane_index]
    plane_rotate = matrix_rotate([in_plane_vector[0], 0, in_plane_vector[2]], [1, 0, 0])
    array_coords = asarray([dot(plane_rotate, x) for x in array_coords])
    return array_coords




fgroup = 'C(=O)ONC'

attached = '[cH]1[cH][cH][cH][cH]c1'

pybel_mol = pybel.readstring('smi', attached + fgroup)
pybel_mol.make3D(forcefield='UFF')

uff = ob.OBForceField_FindForceField('uff')
uff.Setup(pybel_mol.OBMol)
uff.GetAtomTypes(pybel_mol.OBMol)

atom_types = []
coordinates = []

for ob_atom in pybel_mol:
    coordinates.append(ob_atom.coords)
#    atom_types.append(ob_atom.OBAtom.GetData("FFAtomType").GetValue())

rotated_coordinates = realign(coordinates, 11, 10, 8)

bonds = {}

# look at all the bonds separately from the atoms
for bond in ob.OBMolBondIter(pybel_mol.OBMol):
    # These rules are translated from ob/forcefielduff.cpp...
    start_idx = bond.GetBeginAtomIdx()
    end_idx = bond.GetEndAtomIdx()

    start_atom = bond.GetBeginAtom()
    end_atom = bond.GetEndAtom()

    bond_order = bond.GetBondOrder()
    if bond.IsAromatic():
        bond_order = 1.5
    # e.g., in Cp rings, may not be "aromatic" by OB
    # but check for explicit hydrogen counts
    #(e.g., biphenyl inter-ring is not aromatic)
    #FIXME(tdaff): aromatic C from GetType is "Car" is this correct?
    if start_atom.GetType()[-1] == 'R' and end_atom.GetType()[-1] == 'R' and start_atom.ExplicitHydrogenCount() == 1 and end_atom.ExplicitHydrogenCount() == 1:
        bond_order = 1.5
    if bond.IsAmide():
        bond_order = 1.41
    # save the indicies as zero based
    bond_length = bond.GetLength()
    bond_id = tuple(sorted((start_idx-12, end_idx-12)))
    bonds[bond_id] = (bond_length, bond_order)

print bonds

output_text = []

atom_block = []
for ob_atom, coord in zip(pybel_mol, rotated_coordinates):
    element = ATOMIC_NUMBER[ob_atom.atomicnum]
    ff_type = ob_atom.OBAtom.GetData("FFAtomType").GetValue()
    atom_block.append("    {0:4} {1:4} {2[0]:10.6f} {2[1]:10.6f} {2[2]:10.6f}\n".format(element, ff_type, coord))

output_text.append('atoms =\n')
output_text.extend(atom_block[11:])
output_text.extend('orientation = 0.0 1.0 0.0\n')
output_text.append('normal = 0.0 0.0 1.0\n')
output_text.append('carbon_bond = {:.3f}\n'.format(bonds[(-1, 0)][0]))

print pybel_mol.write(format='ascii', opt={'a': 2, 'w': 70})

bonds_block = []
# no bonds < idx 11
for bond in sorted(bonds):
    if not bond[0] < 0 and not bond[1] < 0:
        bonds_block.append("    {0[0]:4} {0[1]:4} {1[1]:5.1f}\n".format(bond, bonds[bond]))

output_text.append('bonds =\n')
output_text.extend(bonds_block[:])


print("".join(output_text))
