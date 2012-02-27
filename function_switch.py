#!/usr/bin/env python

"""
function_switch.py

Alter a known structure with new functional groups ready for fapping.

"""


import copy
import ConfigParser
import re
from math import factorial
from itertools import chain, combinations

import numpy as np
from numpy import array, asarray, dot, cross
from numpy.linalg import norm

from faps import Structure, Atom
from faps import unique
from config import Options


class ModifiableStructure(Structure):
    """
    Derivative of Structure with methods to facilitate function group
    switching. Use as a staging area for new methods to refactor
    into the parent Structure class.

    """

    def gen_normals(self):
        """
        Calculate the normal to the vectors between the first two neighbours.

        """
        self.gen_neighbour_list()

        cell = self.cell.cell
        inv_cell = self.cell.inverse
        cpositions = [atom.ipos(cell, inv_cell) for atom in self.atoms]
        fpositions = [atom.ifpos(inv_cell) for atom in self.atoms]

        # Speed is not so important here but these are derived from the
        # surface ares distance calcualtions so they are explicitly given
        # the fractional positions too
        for at_idx, atom in enumerate(self.atoms):
            conex = [x[1] for x in atom.neighbours[:2]]
            left = min_vect(cpositions[at_idx], fpositions[at_idx],
                            cpositions[conex[0]], fpositions[conex[0]],
                            cell)
            right = min_vect(cpositions[at_idx], fpositions[at_idx],
                             cpositions[conex[1]], fpositions[conex[1]],
                             cell)
            atom.normal = cross(left, right)
        print [atom.normal for atom in self.atoms]


    def gen_connection_table(self):
        """
        Determine connectivity from the neighbour list.

        """
        self.gen_neighbour_list()
        self.attachments = []
        for atom in self.atoms:
            atom.bonds = []
            for neighbour in atom.neighbours:
                if neighbour[0] < 0.5*(atom.vdw_radius + self.atoms[neighbour[1]].vdw_radius):
                    atom.bonds.append(neighbour)
                    break
            if atom.type == 'H':
                neighbour_a = self.atoms[self.atoms[atom.bonds[0][1]].neighbours[0][1]].pos
                neighbour_b = self.atoms[self.atoms[atom.bonds[0][1]].neighbours[1][1]].pos
                pivot = self.atoms[atom.bonds[0][1]].pos
                normal = cross(direction3d(pivot, neighbour_a), direction3d(pivot, neighbour_b))
                self.attachments.append((self.atoms.index(atom), atom.bonds, normal))

    def gen_site_connection_table(self):
        u_atoms = unique(self.atoms, key=lambda x: x.site)
        self.gen_neighbour_list()
        self.attachments = []
        for atom in u_atoms:
            atom.bonds = []
            for neighbour in atom.neighbours:
                if neighbour[0] < 0.5*(atom.vdw_radius + self.atoms[neighbour[1]].vdw_radius):
                    atom.bonds.append(neighbour)
                    break
            if atom.type == 'H':
                self.attachments.append((self.atoms.index(atom), atom.bonds))


class FunctionalGroupLibrary(dict):
    """
    Container for all the available functional groups just subclasses
    the standard dict adding some new methods.

    """
    def from_file(self, library_file_name='functional_groups.lib'):
        """Parse groups from the ConfigParser .ini style file."""
        # just a standard ConfigParser conversion to a dict of
        # FunctionalGroup objects
        library_file = ConfigParser.SafeConfigParser()
        library_file.read(library_file_name)
        for group_name in library_file.sections():
            self[group_name] = FunctionalGroup(library_file.items(group_name))


class FunctionalGroup(object):
    """
    Substitutable functional group. Bunch of atoms with information on how
    to connect to the framework

    """

    def __init__(self, items):
        """Initialize from a list of tuples as the attributes."""
        # These are defaults, best that they are overwritten
        self.atoms = []
        self.orientation = [0, 1, 0]

        # pop the items from a dict giving neater code
        items = dict(items)

        self._parse_atoms(items.pop('atoms'))
        self.orientation = string_to_tuple(items.pop('orientation'), float)
        self.normal = string_to_tuple(items.pop('normals'), float)
        # Arbitrary attributes can be set
        self.__dict__.update(items)

    def _parse_atoms(self, atom_text):
        """Read atom information from the file."""
        self.atoms = []
        for atom in atom_text.splitlines():
            atom = atom.strip().split()
            if not atom:
                continue
            self.atoms.append(Atom(atom[0], [float(x) for x in atom[1:4]]))

    def atoms_attached_to(self, point, bond_length, direction, normal):
        """Return a list of atoms at the specified position."""
        new_atoms = [copy.copy(atom) for atom in self.atoms]
        rotate_matrix = matrix_rotate(self.orientation, direction)
        my_rotated_normal = np.dot(rotate_matrix, self.normal)
        orient_matrix = matrix_rotate(my_rotated_normal, normal)
        for atom in new_atoms:
            atom.pos = np.dot(orient_matrix, np.dot(rotate_matrix, atom.pos))
            atom.pos = (atom.pos + point + bond_length*np.array(direction))
        return new_atoms


def to_xyz(atoms, charges=True):
    """
    Return a list of atoms in xyz compatible format (includes charges)
    """

    if charges:
        xyz = ["%-6s" % atom.type +
               "%9.6f %9.6f %9.6f" % tuple(atom.pos) +
               "%11.6f\n" % atom.charge for atom in atoms]
    else:
        xyz = ["%-6s" % atom.type +
               "%9.6f %9.6f %9.6f" % tuple(atom.pos) for atom in atoms]

    # header needs the atom counts
    xyz = [" %i\n" % len(atoms), "Structure file generated by faps\n"] + xyz

    return xyz


def matrix_rotate(source, target):
    """Create a rotation matrix that will rotate source on to target."""
    # Normalise so there is no scaling in the array
    source = asarray(source)/norm(source)
    target = asarray(target)/norm(target)
    v = cross(source, target)
    c = dot(source, target)
    h = (1 - c)/(np.dot(v, v))
    return array([[c + h*v[0]*v[0], h*v[0]*v[1] - v[2], h*v[0]*v[2] + v[1]],
                  [h*v[0]*v[1] + v[2], c + h*v[1]*v[1], h*v[1]*v[2] - v[0]],
                  [h*v[0]*v[2] - v[1], h*v[1]*v[2] + v[0], c + h*v[2]*v[2]]])


def direction3d(source, target):
    """Return the vector connecting two 3d points."""
    return [target[0] - source[0],
            target[1] - source[1],
            target[2] - source[2]]

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def string_to_tuple(value, dtype=None):
    """Parse a list of items, ignoring whitespace, brackets and commas."""
    value = [x for x in re.split('[\s,\(\)\[\]]*', value) if x]
    if dtype is not None:
        return tuple([dtype(x) for x in value])
    else:
        return tuple(value)


def min_vect(c_coa, f_coa, c_cob, f_cob_in, box):
    """Calculate the closest distance assuming fractional, in-cell coords."""
    f_cob = f_cob_in[:]
    fdx = f_coa[0] - f_cob[0]
    if fdx < -0.5:
        f_cob[0] -= 1
    elif fdx > 0.5:
        f_cob[0] += 1
    fdy = f_coa[1] - f_cob[1]
    if fdy < -0.5:
        f_cob[1] -= 1
    elif fdy > 0.5:
        f_cob[1] += 1
    fdz = f_coa[2] - f_cob[2]
    if fdz < -0.5:
        f_cob[2] -= 1
    elif fdz > 0.5:
        f_cob[2] += 1
    if f_cob == f_cob_in:
        # if nothing has changed, use initial values
        return direction3d(c_coa, c_cob)
    else:
        new_b = [f_cob[0]*box[0][0] + f_cob[1]*box[1][0] + f_cob[2]*box[2][0],
                 f_cob[0]*box[0][1] + f_cob[1]*box[1][1] + f_cob[2]*box[2][1],
                 f_cob[0]*box[0][2] + f_cob[1]*box[1][2] + f_cob[2]*box[2][2]]
        return direction3d(c_coa, new_b)


def main():
    """
    Run the substitution for an input structure.

    """

    job_options = Options()
    job_name = job_options.get('job_name')

    irmof1 = ModifiableStructure("IMOF1")
    irmof1.from_cif("MOF5/277428.cif")

    irmof1.gen_connection_table()
    #irmof1.symmetry_conections()
    irmof1.gen_normals()

    f_groups = FunctionalGroupLibrary()
    f_groups.from_file()
    print(f_groups)

    print len(irmof1.attachments)

    for atom in irmof1.atoms:
        print atom.bonds

    for r in range(len(irmof1.attachments)+1):
        print(factorial(len(irmof1.attachments))/(factorial(r)*factorial((len(irmof1.attachments)-r))))

    for idx, attach_points in enumerate(powerset(irmof1.attachments)):
        attachment = f_groups['NH2']
    #    print attach_points
        new_mof = list(irmof1.atoms)
        for this_point in attach_points:
            attach_id = this_point[0]
            attach_to = this_point[1][0][1]
            attach_at = irmof1.atoms[attach_to].pos
            attach_towards = direction3d(attach_at, irmof1.atoms[attach_id].pos)
            bond_length = 0.5*(irmof1.atoms[attach_to].vdw_radius + attachment.atoms[0].vdw_radius)
            new_mof[attach_id:attach_id+1] = [None]
            #attach_normal = this_point[2]
            attach_normal = irmof1.atoms[attach_to].normal
            new_mof.extend(attachment.atoms_attached_to(attach_at, bond_length, attach_towards, attach_normal))
        new_mof = [an_atom for an_atom in new_mof if an_atom is not None]
        print idx, '\r',
        if not idx % 500:
            with open('output%04i.xyz' % idx, 'wb') as output_file:
                output_file.writelines(to_xyz(new_mof))

if __name__ == '__main__':

    main()
