#!/usr/bin/env python

"""
function_switch.py

Alter a known structure with new functional groups ready for fapping.

"""


import copy
import ConfigParser
import hashlib
import pickle
import random
import re
import sys
import textwrap
import time
from itertools import chain, combinations, product
from logging import debug, info, warn, error
from os import path

import numpy as np
import openbabel as ob
import pybel
from numpy import array, identity, asarray, dot, cross, outer, sin, cos
from numpy.linalg import norm

from faps import Structure, Atom
from faps import vecdist3, subgroup
from config import Options
from elements import UFF_TYPES, CCDC_BOND_ORDERS


class ModifiableStructure(Structure):
    """
    Derivative of Structure with methods to facilitate function group
    switching. Use as a staging area for new methods to refactor
    into the parent Structure class.

    """

    def gen_factional_positions(self):
        """
        Precalculate the fractional positions for all the atoms in the
        current cell. These will be incorrect if the cell chagnes!

        """

        cell = self.cell.cell
        inv_cell = self.cell.inverse
        for atom in self.atoms:
            atom.cellpos = atom.ipos(cell, inv_cell)
            atom.cellfpos = atom.ifpos(inv_cell)

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


    def gen_site_connection_table(self):
        """
        Find connected atoms (bonds) and organise by sites. Symmetry
        equivalent atoms have the same site.

        """

        atoms = self.atoms
        self.gen_neighbour_list()
        self.attachments = {}

        for at_idx, atom in enumerate(self.atoms):
            atom.bonds = []
            neighbours = atom.neighbours

            max_bonds = []
            for neighbour in neighbours:
                n_atom = self.atoms[neighbour[1]]
                if neighbour[0] < (atom.bond_cutoff + n_atom.bond_cutoff):
                    max_bonds.append(neighbour[1])
                elif (neighbour[0] - atom.bond_cutoff) > 2.88:
                    # some non bonding atoms might be closer than bonding ones
                    # but francium has largest cutoff; can stop after 2.88
                    break

            if atom.type == 'H':
                # check if second neighbour is within bridging cutoff
                if neighbours[1][0] <= 1.1*(0.460 + atoms[neighbours[1][1]].bond_cutoff):
                    atom.uff_type = 'H_b'
                    atom.bonds = [neighbours[0][1], neighbours[1][1]]
                else:
                    # assume lone hydrogen
                    atom.uff_type = 'H_'
                    atom.bonds = [neighbours[0][1]]
            elif atom.type == 'B':
                # check if tetrahedral sp3 or sp2
                if neighbours[3][0] <= 1.1*(0.838 + atoms[neighbours[3][1]].bond_cutoff):
                    atom.uff_type = 'B_3'
                    atom.bonds = [x[1] for x in neighbours[:4]]
                else:
                    # planar sp2
                    atom.uff_type = 'B_2'
                    atom.bonds = [x[1] for x in neighbours[:3]]
            elif atom.type == 'C':
                # Four types of Carbon, but C_R and C_2 same params
                if neighbours[3][0] <= 1.1*(0.757 + atoms[neighbours[3][1]].bond_cutoff):
                    atom.uff_type = 'C_3'
                    atom.bonds = [x[1] for x in neighbours[:4]]
                elif neighbours[2][0] <= 1.1*(0.732 + atoms[neighbours[2][1]].bond_cutoff):
                    # planar sp2 and aromatic have the same UFF parameters
                    # TODO(tdaff): make aromatic if necessary
                    atom.uff_type = 'C_2'
                    atom.bonds = [x[1] for x in neighbours[:3]]
                else:
                    # linear sp
                    atom.uff_type = 'C_1'
                    atom.bonds = [x[1] for x in neighbours[:2]]
            elif atom.type == 'N':
                # Careful, N_R and N_2 have different angles
                if neighbours[2][0] <= 1.1*(0.700 + atoms[neighbours[2][1]].bond_cutoff):
                    atom.uff_type = 'N_3'
                    atom.bonds = [x[1] for x in neighbours[:3]]
                elif neighbours[1][0] <= 1.1*(0.699 + atoms[neighbours[1][1]].bond_cutoff):
                    # FIXME(tdaff): make aromatic if necessary
                    atom.uff_type = 'N_2'
                    atom.bonds = [x[1] for x in neighbours[:2]]
                else:
                    # sp
                    atom.uff_type = 'N_1'
                    atom.bonds = [x[1] for x in neighbours[:1]]
            elif atom.type == 'O':
                # UFF has sp3 to sp1 and aromatic oxygen
                # FIXME(tdaff): might need aromatic check
                if neighbours[1][0] <= 1.1*(0.680 + atoms[neighbours[1][1]].bond_cutoff):
                    atom.uff_type = 'O_3'
                    atom.bonds = [x[1] for x in neighbours[:2]]
                else:
                    # terminal oxygen?
                    atom.uff_type = 'O_2'
                    atom.bonds = [x[1] for x in neighbours[:1]]
            elif atom.type == 'P':
                atom.bonds = max_bonds
                if len(max_bonds) > 3:
                    atom.uff_type = "P_3+5"
                else:
                    atom.uff_type = "P_3+3"
            elif atom.type == 'S':
                warn("Sulphr not checked")
                atom.uff_type = "S_3"
                atom.bonds = [x[1] for x in neighbours[:4]]

            else:
                atom.uff_type = UFF_TYPES[atom.type]
                atom.bonds = max_bonds

            if atom.type == 'H':
                if atom.site in self.attachments:
                    self.attachments[atom.site].append((self.atoms.index(atom), atom.bonds))
                else:
                    self.attachments[atom.site] = [(self.atoms.index(atom), atom.bonds)]

    def gen_babel_uff_properties(self):
        """
        Process a supercell with babel to calculate UFF atom types and
        bond orders.
        """
        # Pass as free form fractional
        # make a 2x2x2 supercell with the original atom in the centre
        cell = self.cell
        super_cell = (cell.a*2, cell.b*2, cell.c*2,
                      cell.alpha, cell.beta, cell.gamma)
        as_fffract = ['generated fractionals\n', '%f %f %f %f %f %f\n' % super_cell]
        for x_image in [0, 1]:
            for y_image in [0, 1]:
                for z_image in [0, 1]:
                    for atom in self.atoms:
                        ifpos = atom.cellfpos
                        if ifpos[0] < 0.5:
                            new_xpos = (ifpos[0] + x_image)/2.0
                        else:
                            new_xpos = (ifpos[0] - x_image)/2.0
                        if ifpos[1] < 0.5:
                            new_ypos = (ifpos[1] + y_image)/2.0
                        else:
                            new_ypos = (ifpos[1] - y_image)/2.0
                        if ifpos[2] < 0.5:
                            new_zpos = (ifpos[2] + z_image)/2.0
                        else:
                            new_zpos = (ifpos[2] - z_image)/2.0
                        ifpos = [new_xpos, new_ypos, new_zpos]
                        atom_line = ("%s " % atom.type +
                                     "%f %f %f\n" % tuple(ifpos))
                        as_fffract.append(atom_line)
        pybel_string = ''.join(as_fffract)
        pybel_mol = pybel.readstring('fract', pybel_string)
        # need to tell the typing system to ignore all atoms in the setup
        # or it will silently crash with memory issues
        constraint = ob.OBFFConstraints()
        for at_idx in range(pybel_mol.OBMol.NumAtoms()):
            constraint.AddIgnore(at_idx)
        uff = ob.OBForceField_FindForceField('uff')
        uff.Setup(pybel_mol.OBMol, constraint)
        uff.GetAtomTypes(pybel_mol.OBMol)
        for atom, ob_atom in zip(self.atoms, pybel_mol):
            atom.uff_type = ob_atom.OBAtom.GetData("FFAtomType").GetValue()

        bonds = {}
        max_idx = self.natoms
        # look at all the bonds separately from the atoms
        for bond in ob.OBMolBondIter(pybel_mol.OBMol):
            # These rules are translated from ob/forcefielduff.cpp...
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if start_idx > max_idx and end_idx > max_idx:
                continue
            if end_idx > max_idx:
                end_idx = end_idx % max_idx
            if start_idx > max_idx:
                start_idx = start_idx % max_idx

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
            bonds[tuple(sorted((start_idx-1, end_idx-1)))] = bond_order

        self.bonds = bonds


def to_cif(atoms, cell, bonds, name):
    """Return a CIF file with bonding and atom types."""

    inv_cell = cell.inverse

    type_count = {}

    atom_part = []
    for idx, atom in enumerate(atoms):
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

    bond_part = []
    for bond, order in bonds.items():
        try:
            bond_part.append("%-5s %-5s %-5s\n" %
                             (atoms[bond[0]].site, atoms[bond[1]].site,
                              CCDC_BOND_ORDERS[order]))
        except AttributeError:
            # one of the atoms is None so skip
            debug("cif NoneType atom")

    cif_file = [
        "data_%s\n" % name.replace(' ', '_'),
        "%-33s %s\n" % ("_audit_creation_date", time.strftime('%Y-%m-%dT%H:%M:%S%z')),
        "%-33s %s\n" % ("_audit_creation_method", "LUBE"),
        "%-33s %s\n" % ("_symmetry_space_group_name_H-M", "P1"),
        "%-33s %s\n" % ("_symmetry_Int_Tables_number", "1"),
        "%-33s %s\n" % ("_space_group_crystal_system", cell.crystal_system),
        "%-33s %s\n" % ("_cell_length_a", cell.a),
        "%-33s %s\n" % ("_cell_length_b", cell.b),
        "%-33s %s\n" % ("_cell_length_c", cell.c),
        "%-33s %s\n" % ("_cell_angle_alpha", cell.alpha),
        "%-33s %s\n" % ("_cell_angle_beta", cell.beta),
        "%-33s %s\n" % ("_cell_angle_gamma", cell.gamma),
        # start of atom loops
        "\nloop_\n",
        "_atom_site_label\n",
        "_atom_site_type_symbol\n",
        "_atom_type_description\n",
        "_atom_site_fract_x\n",
        "_atom_site_fract_y\n",
        "_atom_site_fract_z\n"] + atom_part + [
        # bonding loop
        "\nloop_\n",
        "_geom_bond_atom_site_label_1\n",
        "_geom_bond_atom_site_label_2\n",
#        "_geom_bond_distance\n",
        "_ccdc_geom_bond_type\n"] + bond_part

    return cif_file


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
        debug("Reading groups from %s" % library_file_name)
        library_file.read(library_file_name)
        for group_name in library_file.sections():
            self[group_name] = FunctionalGroup(library_file.items(group_name))

    @property
    def group_list(self):
        """Compile a soretd list of all the available groups."""
        return sorted(list(self))


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
        self.normal = string_to_tuple(items.pop('normal'), float)
        self.bond_length = float(items.pop('carbon_bond'))
        self.bonds = dict(((int(x), int(y)), float(z)) for (x, y, z) in subgroup(items.pop('bonds').split(), width=3))
        self.idx = 0
        self.connection_point = 0  # always connect to the first atom
        # Arbitrary attributes can be set
        self.__dict__.update(items)
        self._gen_neighbours()

    def _parse_atoms(self, atom_text):
        """Read atom information from the file."""
        self.atoms = []
        for atom in atom_text.splitlines():
            atom = atom.strip().split()
            if not atom:
                continue
            new_atom = Atom(atom[0], [float(x) for x in atom[2:5]])
            new_atom.uff_type = atom[1]
            new_atom.site = label_atom(new_atom.element)
            self.atoms.append(new_atom)

    def _gen_neighbours(self):
        """Update atoms with neighbouring atoms."""
        # Iterate over all pairs
        # Assume non periodic and geometries are good
        # dummy atom at the tether point
        self.atoms.append(Atom('C', [-self.bond_length*x for x in self.orientation]))
        for atom in self.atoms:
            # distance matrix for all neighbours except self
            neighbours = []
            for ot_idx, other in enumerate(self.atoms):
                if atom is other:
                    continue
                length = vecdist3(atom.pos, other.pos)
                neighbours.append((length, ot_idx))
            # these are expected to be is distance order
            atom.neighbours = sorted(neighbours)
        # remove the dummy atom
        self.atoms.pop()
        # first atom should be bonded to the tether, flag this as '-1'
        #self.atoms[0].bonds[self.atoms[0].bonds.index(self.natoms)] = -1

    def atoms_attached_to(self, point, direction, normal, attach_point, start_index, bond_length=None):
        """Return a list of atoms at the specified position."""
        if bond_length is None:
            bond_length = self.bond_length
        new_atoms = [copy.copy(atom) for atom in self.atoms]
        rotate_matrix = matrix_rotate(self.orientation, direction)
        my_rotated_normal = np.dot(rotate_matrix, self.normal)
        orient_matrix = matrix_rotate(my_rotated_normal, normal)
        for index, atom in enumerate(new_atoms):
            atom.pos = np.dot(orient_matrix, np.dot(rotate_matrix, atom.pos))
            atom.pos = (atom.pos + point + bond_length*np.array(direction))
            atom.idx = start_index + index
        new_bonds = {}
        for bond_pair, bond_order in self.bonds.iteritems():
            new_bond = (bond_pair[0] + start_index, bond_pair[1] + start_index)
            new_bonds[new_bond] = bond_order
        # bond to structure
        new_bonds[(attach_point, self.connection_point + start_index)] = 1

        return new_atoms, new_bonds

    @property
    def natoms(self):
        """The number of atoms in the functional group."""
        return len(self.atoms)


def to_pdb(atoms, cell, name=None):
    """
    Return a pdb compatible file representation.

    """
    pdb = [
        "TITLE  faps functional group switching\n",
        "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           0\n" % cell.params]
    site_number = 1
    for at_idx, atom in enumerate(atoms):
        if atom.site is None:
            atom.site = "%s%i" % (atom.type, site_number)
            site_number += 1
        pdb.append("ATOM  %5i %-4s UNK A   1    " % (at_idx+1, atom.site) +
                   "%8.3f%8.3f%8.3f" % tuple(atom.pos) +
                   #"%6.2f  0.00%12s  \n" % (atom.charge, atom.type))
                   # charge -> occupancy in vesta
                   "%6.2f  0.00%12s  \n" % (1.0, atom.type))

    if name is not None:
        pdb.extend(["REMARK".ljust(70), "\n",
                    textwrap.fill(name, initial_indent="REMARK     ",
                                  subsequent_indent="REMARK     ",
                                  break_long_words=True), "\n",
                    "REMARK".ljust(70), "\n"])

    return pdb


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


def rotation_about_angle(axis_in, angle):
    """
    Create a rotation matrix corresponding to the rotation around a general
    axis by a specified angle.

    """
    axis = axis_in/norm(axis_in)

    aprod = outer(axis, axis)
    skew = array([[0, axis[2], -axis[1]], [-axis[2], 0, axis[0]],
        [axis[1], -axis[0], 0]])

    # R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)
    return aprod+cos(angle)*(identity(3)-aprod)+sin(angle)*skew


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


def test_collision(test_atom, atoms, cell, overlap=1.3):
    """
    Does atom intersect with any others?

    """
    pos = test_atom.ipos(cell.cell, cell.inverse)
    ipos = test_atom.ifpos(cell.inverse)
    for atom in atoms:
        if atom is None:
            continue
        dist = min_vect(pos, ipos, atom.ipos(cell.cell, cell.inverse), atom.ifpos(cell.inverse), cell.cell)
        dist = dot(dist, dist)
        if dist < overlap:
            return False
    return True


def uff_bonds(atom, structure, tolerance=1.01):
    atom.bonds = []
    max_bonds = []
    neighbours = atom.neighbours
    atoms = structure.atoms

    for neighbour in neighbours:
        n_atom = atoms[neighbour[1]]
        if neighbour[0] < tolerance*(atom.bond_cutoff + n_atom.bond_cutoff):
            max_bonds.append(neighbour[1])
        elif (neighbour[0] - atom.bond_cutoff) > 2.88:
            # some non bonding atoms might be closer than bonding ones
            # but francium has largest cutoff; can stop after 2.88
            break

    if atom.type == 'H':
        # check if second neighbour is within bridging cutoff
        if len(max_bonds) > 1:
            atom.uff_type = 'H_b'
            atom.bonds = max_bonds[:2]
        else:
            # assume lone hydrogen
            atom.uff_type = 'H_'
            atom.bonds = max_bonds
    elif atom.type == 'B':
        # check if tetrahedral sp3 or sp2
        if len(max_bonds) > 3:
            atom.uff_type = 'B_3'
            atom.bonds = max_bonds[:4]
        else:
            # planar sp2
            atom.uff_type = 'B_2'
            atom.bonds = max_bonds[:3]
    elif atom.type == 'C':
        # Four types of Carbon, but C_R and C_2 same params
        if len(max_bonds) > 3:
            atom.uff_type = 'C_3'
            atom.bonds = max_bonds[:4]
        elif len(max_bonds) > 2:
            # planar sp2 and aromatic have the same UFF parameters
            # TODO(tdaff): make aromatic if necessary
            atom.uff_type = 'C_2'
            atom.bonds = max_bonds[:3]
        else:
            # linear sp
            atom.uff_type = 'C_1'
            atom.bonds = max_bonds[:2]
    elif atom.type == 'N':
        # Careful, N_R and N_2 have different angles
        if len(max_bonds) > 2:
            atom.uff_type = 'N_3'
            atom.bonds = max_bonds[:3]
        elif len(max_bonds) > 1:
            # FIXME(tdaff): make aromatic if necessary
            atom.uff_type = 'N_2'
            atom.bonds = max_bonds[:2]
        else:
            # sp
            atom.uff_type = 'N_1'
            atom.bonds = max_bonds[:1]
    elif atom.type == 'O':
        # UFF has sp3 to sp1 and aromatic oxygen
        # FIXME(tdaff): might need aromatic check
        if len(max_bonds) > 1:
            atom.uff_type = 'O_3'
            atom.bonds = max_bonds[:2]
        else:
            # terminal oxygen?
            atom.uff_type = 'O_2'
            atom.bonds = max_bonds[:1]
    elif atom.type == 'P':
        if len(max_bonds) > 3:
            atom.uff_type = "P_3+5"
            atom.bonds = max_bonds
        else:
            atom.uff_type = "P_3+3"
            atom.bonds = max_bonds[:3]
    elif atom.type == 'S':
        warn("Sulphur not checked")
        atom.uff_type = "S_3"
        atom.bonds = max_bonds[:4]
    else:
        atom.uff_type = UFF_TYPES[atom.type]
        atom.bonds = max_bonds



def all_combinations_replace(structure, groups, rotations=12, replace_only=None):
    """
    Replace every functional point with every combination of functional groups.

    """
    if replace_only is not None:
        local_attachments = [att_id for att_id in structure.attachments if att_id in replace_only]
        debug("Replacing only: %s" % list(local_attachments))
    else:
        local_attachments = structure.attachments
        debug("Replacing all sites: %s" % local_attachments.keys())
    sites = powerset(sorted(local_attachments))
    for site_set in sites:
        for group_set in product(groups, repeat=len(site_set)):
            replace_list = zip(group_set, site_set)
            site_replace(structure, groups, replace_list, rotations=12)

def site_replace(structure, groups, replace_list, rotations=12):
    """
    Replace atoms at site_set with corresponding items from group_set.

    Will write files on success, and return 1 for failed attempt.

    """
    rotation_angle = 2*np.pi/rotations
    new_mof_name = []
    new_mof_friendly_name = []
    # copy the atoms and bonds so we don't alter the original structure
    new_mof = list(structure.atoms)
    new_mof_bonds = dict(structure.bonds)
    for this_group, this_site in replace_list:
        attachment = groups[this_group]
        new_mof_name.append("%s@%s" % (this_group, this_site))
        new_mof_friendly_name.append("%s@%s" % (attachment.name, this_site))
        for this_point in structure.attachments[this_site]:
            attach_id = this_point[0]
            attach_to = this_point[1][0]
            attach_at = structure.atoms[attach_to].pos
            attach_towards = direction3d(attach_at, structure.atoms[attach_id].pos)
            attach_normal = structure.atoms[attach_to].normal
            new_mof[attach_id:attach_id+1] = [None]
            start_idx = len(new_mof)
            for _trial_rotation in range(rotations):
                incoming_group, incoming_bonds = attachment.atoms_attached_to(attach_at, attach_towards, attach_normal, attach_to, start_idx)
                for atom in incoming_group:
                    if not test_collision(atom, new_mof, structure.cell):
                        debug("Rotating group")
                        attach_normal = dot(rotation_about_angle(attach_towards, rotation_angle), attach_normal)
                        # Don't need to test more atoms
                        break
                else:
                    # Fits, so add and move on
                    new_mof.extend(incoming_group)
                    new_mof_bonds.update(incoming_bonds)
                    break
            else:
                # Did not attach
                error("Failed: %s from %s" % (this_group, group_set))
                return 1
    new_mof_name = ".".join(new_mof_name)
    new_mof_friendly_name = ".".join(new_mof_friendly_name)
    info("Generated (%i): [%s]" % (count(), new_mof_friendly_name))
    job_name = structure.name
    with open('%s_func_%s.cif' % (job_name, new_mof_name), 'w') as output_file:
        output_file.writelines(to_cif(new_mof, structure.cell, new_mof_bonds, new_mof_name))
    new_mof = [an_atom for an_atom in new_mof if an_atom is not None]
    with open('%s_func_%s.pdb' % (job_name, new_mof_name), 'w') as output_file:
        output_file.writelines(to_pdb(new_mof, structure.cell, name=new_mof_name))


def random_replace(structure, groups, num_groups=None, custom=None, rotations=36):
    """
    Replace a random number of sites.

    """
    nsites = sum(len(x) for x in structure.attachments.values())
    if custom is not None:
        debug("Processing custom string: %s" % custom)
        func_repr = custom.strip('{}').split(".")
        if len(func_repr) != nsites:
            error("Expected %s sites; got %s" % (nsites, len(func_repr)))
    else:
        if num_groups is None:
            num_groups = random.randint(1, nsites)
            debug("Randomly replacing %i sites" % num_groups)
        elif num_groups > nsites:
            warn("Too many sites requested; changing all %i" % nsites)
            num_groups = nsites
        #TODO(tdaff): selected groups only
        func_repr = [random.choice(groups.keys()) for _counter in range(num_groups)]
        # Pad to the correct length
        func_repr.extend([""]*(nsites - num_groups))
        # Randomise
        random.shuffle(func_repr)
    # Unique-ish
    unique_name = hashlib.md5(str(func_repr)).hexdigest()
    new_mof_name = []
    new_mof = list(structure.atoms)
    new_mof_bonds = dict(structure.bonds)
    for this_point, this_group in zip(chain(*[structure.attachments[x] for x in sorted(structure.attachments)]), func_repr):
        if this_group == "":
            new_mof_name.append("")
            continue
        else:
            new_mof_name.append(this_group)
        attachment = groups[this_group]
        attach_id = this_point[0]
        attach_to = this_point[1][0]
        attach_at = structure.atoms[attach_to].pos
        attach_towards = direction3d(attach_at, structure.atoms[attach_id].pos)
        attach_normal = structure.atoms[attach_to].normal
        #extracted_atoms = new_mof[attach_id:attach_id+1]
        new_mof[attach_id:attach_id+1] = [None]
        start_idx = len(new_mof)
        for trial_rotation in range(rotations):
            incoming_group, incoming_bonds = attachment.atoms_attached_to(attach_at, attach_towards, attach_normal, attach_to, start_idx)
            for atom in incoming_group:
                if not test_collision(atom, new_mof, structure.cell):
                    debug("Randomly rotating group")
                    attach_normal = dot(rotation_about_angle(attach_towards, random.random()*np.pi*2), attach_normal)
                    break
            else:
                # Fits, so add and move on
                new_mof.extend(incoming_group)
                new_mof_bonds.update(incoming_bonds)
                break
        else:
            # this_point not valid
            error("Failed to generate: %s" % ".".join([x or "" for x in func_repr]))
            warn("Stopped after: %s" % ".".join(new_mof_name))

            return False

    new_mof_name = "{" + ".".join(new_mof_name) + "}"
    info("Generated (%i): %s" %  (count(), new_mof_name))
    info("With unique name: %s" % unique_name)
    with open('random-%s.cif' % unique_name, 'wb') as output_file:
        output_file.writelines(to_cif(new_mof, structure.cell, new_mof_bonds, new_mof_name))
    new_mof = [an_atom for an_atom in new_mof if an_atom is not None]
    with open('random-%s.pdb' % unique_name, 'wb') as output_file:
        output_file.writelines(to_pdb(new_mof, structure.cell, name=new_mof_name))

    # completed sucessfully
    return True


def label_atom(element=None, site=None):
    """Produce unique atom labels."""
    if not hasattr(label_atom, 'seen'):
        label_atom.seen = set()
    if not hasattr(label_atom, 'index'):
        label_atom.index = 1
    if element is not None:
        label = "%s%i" % (element, label_atom.index)
        while label in label_atom.seen:
            label_atom.index += 1
            label = "%s%i" % (element, label_atom.index)
        label_atom.seen.add(label)
        return label
    elif site is not None:
        label_atom.seen.add(site)


def count(reset=False):
    """Return the next itneger from a global state."""
    if not hasattr(count, 'idx') or reset is True:
        count.idx = 0
    elif reset:
        count.idx = reset
    else:
        count.idx += 1
    return count.idx


def main():
    """
    Run the substitution for an input structure.

    """

    job_options = Options()
    job_name = job_options.get('job_name')

    # Load an existing pickled structure or generate a new one
    pickle_file = "__%s.lube_structure" % job_name
    if path.exists(pickle_file):
        info("Existing structure found: %s; loading..." % pickle_file)
        with open(pickle_file, 'rb') as load_structure:
            input_structure = pickle.load(load_structure)
    else:
        info("Initialising a new structure. This may take some time.")
        input_structure = ModifiableStructure(job_name)
        input_structure.from_file(job_name,
                                  job_options.get('initial_structure_format'),
                                  job_options)

        input_structure.gen_site_connection_table()
        input_structure.gen_normals()

        # Ensure that atoms in the structure are properly typed
        input_structure.gen_factional_positions()
        input_structure.gen_babel_uff_properties()

        # Cache the results
        info("Dumping cache of structure connectivity.")
        with open(pickle_file, 'wb') as save_structure:
            pickle.dump(input_structure, save_structure)

    # label_atom has a global state that
    for atom in input_structure.atoms:
        label_atom(site=atom.site)

    f_groups = FunctionalGroupLibrary()
    f_groups.from_file()
    debug("Groups in library: %s" % str(f_groups.group_list))

    if job_options.getbool('daemon'):

        # We need to process inputs as the come in; regex out the strings
        info("Waiting for user input....")

        line = raw_input('fswitch >>> ')
        while line:
        #for line in sys.stdin.readlines():
            if 'help' in line.lower():
                info("Fapswitch daemon mode:")
                info("Random strings: dot separated between curly braces")
                info("e.g. {.Me...F..Cl.Cl.Cl..}")
                info("Site replacements: dot separated between square brackets")
                info("e.g. [Me@H7.COOH@H8.F@H12]")
                info("Multiple structures can be input on a single line")
                info("Blank line exits interative mode")

            # randoms are in braces {}, no spaces
            randoms = re.findall('\{(.*?)\}', line)
            debug("Random strings: %s" % str(randoms))
            for random_string in randoms:
                random_replace(input_structure, f_groups, custom=random_string)

            # sites need unpairing into two zippable lists
            site_strings = re.findall('\[(.*?)\]', line)
            debug("Site replacement strings: %s" % str(site_strings))
            for site_string in site_strings:
                site_list = [x.split('@') for x in site_string.split('.') if x]
                debug(str(site_list))
                site_replace(input_structure, f_groups, site_list)

            # Next line of input
            line = raw_input('fswitch >>> ')


    # Will use selected groups if specified, otherwise use all
    try:
        replace_only = job_options.gettuple('lube_replace_only')
    except AttributeError:
        replace_only = None
    if replace_only == ():
        replace_only = None

    if job_options.getbool('lube_replace_all'):
        all_combinations_replace(input_structure, f_groups, replace_only=replace_only)

    custom_strings = job_options.gettuple('lube_custom_strings')
    for custom_string in custom_strings:
        random_replace(input_structure, f_groups, custom=custom_string)

    random_count = job_options.getint('lube_random_structure_count')
    successful_randoms = 0
    while successful_randoms < random_count:
        #function returns true if structure is generated
        if random_replace(input_structure, f_groups):
            successful_randoms += 1


if __name__ == '__main__':

    main()
