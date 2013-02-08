# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

import numpy as np
from numpy import *
from numpy.linalg import *

# <codecell>

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

# <codecell>

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

# <codecell>

def xyz_to_array(text):
    names = [x.split()[0] for x in text.splitlines()]
    coords = asarray([[float(y) for y in x.split()[1:4]] for x in text.splitlines()])
    return names, coords

def array_to_xyz(names, in_coords):
    for name, pos in zip(names, in_coords):
        print("%-8s" % name + "%12.6f %12.6f %12.6f" % tuple(pos))

# <markdowncell>

# Ethyl/Propyl
# ------------

# <codecell>

f_names, f_group = xyz_to_array("""C    0.00000000  0.00000000  0.00000000
C   -1.43577200 -0.45451500  0.22172500
C   -2.26467000 -0.74801700 -0.86619100
C   -3.58036000 -1.17226100 -0.68018600
C   -4.09381400 -1.31155200  0.60683100
C   -3.28018800 -1.02378800  1.70235600
C   -1.96707900 -0.60039000  1.51017400
H   -1.35148300 -0.38146700  2.37624700
H   -3.66943300 -1.12893700  2.70933500
H   -5.11614200 -1.64023100  0.75622400
H   -4.20234000 -1.39232300 -1.54109900
H   -1.87481700 -0.64163100 -1.87397500
C    1.00314700 -1.05522200  0.52241500
C    2.46545500 -0.80655200  0.13186900
C    3.40006500 -1.92810600  0.59672000
H    3.37467900 -2.03780200  1.68555200
H    3.11153900 -2.89010800  0.16175100
H    4.43589500 -1.72966900  0.30753900
H    2.53464800 -0.70359000 -0.95833900
H    2.81107400  0.14386000  0.55157100
H    0.92491900 -1.11714500  1.61541500
H    0.69838100 -2.03712200  0.14173500
H    0.14363800  0.07299200 -1.08516100
C    0.24915600  1.39906200  0.59402200
H    0.17767700  1.38577700  1.68588500
H    1.24165300  1.77202300  0.33007800
H   -0.48748100  2.11528500  0.22142000""")

# <codecell>

new_group = realign(f_group, 12, 0, 13)
print(array_to_xyz(f_names, new_group))

# <markdowncell>

# Methoxy
# -------

# <codecell>

f_names, f_group = xyz_to_array("""C    0.00000000  0.00000000  0.00000000
H   -0.04142900  0.63046100  0.89589700
H   -0.04147800  0.63077200 -0.89567600
H    0.92016900 -0.59054100 -0.00013400
O   -0.99341100 -0.99449600 -0.00011500
C   -2.31061300 -0.57885600 -0.00007200
C   -2.74118300  0.75126300 -0.00006900
C   -4.10646900  1.02085400 -0.00006100
C   -5.03170100 -0.01484600 -0.00004000
C   -4.59551700 -1.33802200 -0.00002800
C   -3.24015200 -1.63205800 -0.00004800
H   -2.89805500 -2.67325800 -0.00004500
H   -5.32724700 -2.15263400  0.00000100
H   -6.10343400  0.20699700 -0.00003100
H   -4.44982000  2.06062600 -0.00007900
H   -2.02821700  1.58310800 -0.00006000""")

# <codecell>

new_group = realign(f_group, 4, 5, 3)
array_to_xyz(f_names, new_group)

# <codecell>

f_names, f_group = xyz_to_array("""C    0.00000000  0.00000000  0.00000000
H    0.12208000  0.95256900 -0.53292700
H   -0.47971700  0.22173300  0.96175700
H   -0.70071100 -0.61284000 -0.58271500
C    1.32298300 -0.70238500  0.19879700
H    1.17633000 -1.63935300  0.77132900
H    2.00016400 -0.07936600  0.81680800
C    1.97971500 -1.00305400 -1.14724600
H    1.35404400 -1.66294200 -1.77284300
H    2.17474000 -0.06401400 -1.70974000
O    3.15641800 -1.78017000 -1.00707600
C    4.32496600 -1.10104200 -0.71554600
C    4.51087000  0.27124000 -0.91440300
C    5.74942400  0.83686200 -0.63738200
C    6.79392800  0.04697100 -0.17058400
C    6.60385900 -1.31902000  0.02126800
C    5.37391800 -1.90507500 -0.24590900
H    5.22577100 -2.98042900 -0.09660400
H    7.42996600 -1.93770400  0.38713900
H    7.76710700  0.49870400  0.04569800
H    5.90014100  1.91039600 -0.79061400
H    3.67926100  0.88515600 -1.29021500""")
new_group = realign(f_group, 10, 11, 7)
array_to_xyz(f_names, new_group)

# <codecell>

f_names, f_group = xyz_to_array("""C    0.00000000  0.00000000  0.00000000
H   -0.15542300 -0.63239300  0.88405600
H   -0.78293500  0.76881800  0.00009100
H   -0.15517300 -0.63192800 -0.88443300
C    1.36987600  0.65253400  0.00037000
H    1.51054100  1.29071800 -0.89346100
H    1.51045200  1.28995200  0.89475700
O    2.32009400 -0.40740200  0.00002300
C    3.65212900 -0.05768800  0.00006800
C    4.52475100 -1.15974100  0.00014100
C    5.89387500 -0.94033800  0.00015700
C    6.40191900  0.35695000  0.00011200
C    5.53423900  1.44105900  0.00006900
C    4.15611100  1.24697800  0.00004800
H    3.49326300  2.11918700  0.00003100
H    5.93384500  2.46057200  0.00005400
H    7.48410200  0.52020000  0.00013000
H    6.57969100 -1.79397800  0.00020700
H    4.12701400 -2.18107000  0.00018300""")
new_group = realign(f_group, 7, 8, 4)
array_to_xyz(f_names, new_group)

# <codecell>

f_names, f_group = xyz_to_array("""C    0.00000000  0.00000000  0.00000000
N    1.11803515  0.05136289 -1.02847424
C    1.56813756 -1.22193563 -1.72598318
C    0.94217033 -2.46336995 -1.41581208
C    1.36399328 -3.65666845 -2.06949726
C    2.41178345 -3.60853264 -3.03335355
C    3.03775068 -2.36709832 -3.34352465
C    2.61592773 -1.17379982 -2.68983946
H    3.09490582 -0.22387815 -2.92717642
H    3.83949957 -2.33026573 -4.08104894
H    2.73455426 -4.52162171 -3.53354088
H    0.88501520 -4.60659011 -1.83216030
H    0.14042144 -2.50020254 -0.67828779
H    1.58822464  0.98385479 -1.26145639
H   -0.18222116  1.00181706  0.38891905
H    0.28865914 -0.66029063  0.81779715
H   -0.90818687 -0.37835903 -0.46919191""")
new_group = realign(f_group, 1, 2, 0)
array_to_xyz(f_names, new_group)

# <codecell>

f_names, f_group = xyz_to_array("""S    0.00000000  0.00000000  0.00000000
C   -1.76798600  0.01178600  0.08395600
C   -2.48551600 -1.18650200  0.14425100
C   -3.87398600 -1.18175500  0.20826900
C   -4.56710000  0.02390700  0.21351400
C   -3.86573200  1.22377900  0.15288100
C   -2.47772100  1.21632500  0.08975500
H   -1.93297300  2.16923400  0.03692100
H   -4.40646100  2.17638900  0.15459100
H   -5.66133700  0.02894200  0.26435200
H   -4.42075200 -2.12983100  0.25441500
H   -1.94656200 -2.14411100  0.13502600
O    0.36545400 -0.08588700  1.64420000
H    1.04382100  0.49998900  1.95027800
O    0.48928400 -1.21285000 -0.58712600
O    0.52717600  1.27436900 -0.41938600""")
new_group = realign(f_group, 0, 1, 12)
array_to_xyz(f_names, new_group)

# <codecell>

import openbabel as ob
import pybel

fgroup = 'C(=O)O'

attached = 'c1cccc1'

print pybel._forcefields

pybel_mol = pybel.readstring('smi', attached)
pybel_mol.make3D(forcefield='UFF')

#uff = ob.OBForceField_FindForceField('uff')
#uff.Setup(pybel_mol.OBMol)
#uff.GetAtomTypes(pybel_mol.OBMol)

atom_types = []

for ob_atom in pybel_mol:
    atom_types.append(ob_atom.OBAtom.GetData("FFAtomType").GetValue())

print atom_types

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
    bond_id = tuple(sorted((start_idx-1, end_idx-1)))
    bonds[bond_id] = (bond_length, bond_order)

print bonds
