#!/usr/bin/env python

"""
function_switch.py

Alter a known structure with new functional groups ready for fapping.

"""


import copy
import ConfigParser
import hashlib
import random
import re
import textwrap
from logging import debug, info, warn, error
from math import factorial
from itertools import chain, combinations, product

import numpy as np
from numpy import array, identity, asarray, dot, cross, outer, sin, cos
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
                else:
                    break
            if atom.type == 'H':
                neighbour_a = self.atoms[self.atoms[atom.bonds[0][1]].neighbours[0][1]].pos
                neighbour_b = self.atoms[self.atoms[atom.bonds[0][1]].neighbours[1][1]].pos
                pivot = self.atoms[atom.bonds[0][1]].pos
                normal = cross(direction3d(pivot, neighbour_a), direction3d(pivot, neighbour_b))
                self.attachments.append((self.atoms.index(atom), atom.bonds, normal))

    def gen_site_connection_table(self):
        """
        Find connected atoms (bonds) and organise by sites. Symmetry
        equivalent atoms have the same site.

        """

        self.gen_neighbour_list()
        self.attachments = {}
        for at_idx, atom in enumerate(self.atoms):
            atom.bonds = []
            for neighbour in atom.neighbours:
                n_atom = self.atoms[neighbour[1]]
                if neighbour[0] < 0.5*(atom.vdw_radius + n_atom.vdw_radius):
                    atom.bonds.append((n_atom.type, neighbour[1]))
                else:
                    break
            if atom.type == 'H':
                if atom.site in self.attachments:
                    self.attachments[atom.site].append((self.atoms.index(atom), atom.bonds))
                else:
                    self.attachments[atom.site] = [(self.atoms.index(atom), atom.bonds)]



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
        self.normal = string_to_tuple(items.pop('normal'), float)
        self.bond_length = float(items.pop('carbon_bond'))
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

    def atoms_attached_to(self, point, direction, normal, bond_length=None):
        """Return a list of atoms at the specified position."""
        if bond_length is None:
            bond_length = self.bond_length
        new_atoms = [copy.copy(atom) for atom in self.atoms]
        rotate_matrix = matrix_rotate(self.orientation, direction)
        my_rotated_normal = np.dot(rotate_matrix, self.normal)
        orient_matrix = matrix_rotate(my_rotated_normal, normal)
        for atom in new_atoms:
            atom.pos = np.dot(orient_matrix, np.dot(rotate_matrix, atom.pos))
            atom.pos = (atom.pos + point + bond_length*np.array(direction))
        return new_atoms


def to_xyz(atoms, charges=True, name=None):
    """
    Return a list of atoms in xyz compatible format (includes charges)
    """

    if charges:
        xyz = ["%-6s" % atom.type +
               "%11.6f %11.6f %11.6f" % tuple(atom.pos) +
               "%11.6f\n" % atom.charge for atom in atoms]
    else:
        xyz = ["%-6s" % atom.type +
               "%11.6f %11.6f %11.6f" % tuple(atom.pos) for atom in atoms]

    # header needs the atom counts
    if name is None:
        name = "Structure file generated by faps"
    xyz = [" %i\n" % len(atoms), "%s\n" % name] + xyz

    return xyz


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
                   "%6.2f  0.00%12s  \n" % (atom.charge, atom.type))

    if name is not None:
        pdb.extend(["REMARK".ljust(70), "\n",
                    textwrap.fill(name, initial_indent="REMARK     ",
                                  subsequent_indent="REMARK     ",
                                  break_long_words=True), "\n",
                    "REMARK".ljust(70), "\n",])

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
    skew = array([[ 0, axis[2], -axis[1]], [-axis[2], 0, axis[0]],
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

def all_combinations_replace(structure, groups, rotations=12):
    """
    Replace every functional point with every combination of functional groups.

    """
    sites = powerset(sorted(structure.attachments))
    structure.gen_factional_positions()
    idx = 0
    rotation_angle = 2*np.pi/rotations
    for site_set in sites:
        for group_set in product(groups, repeat=len(site_set)):
            idx += 1
            new_mof_name = []
            new_mof = list(structure.atoms)
            did_not_attach = False
            for this_site, this_group in zip(site_set, group_set):
                new_mof_name.append("<%s-%s>" % (this_site, groups[this_group].name))
                attachment = groups[this_group]
                for this_point in structure.attachments[this_site]:
                    attach_id = this_point[0]
                    attach_to = this_point[1][0][1]
                    attach_at = structure.atoms[attach_to].pos
                    attach_towards = direction3d(attach_at, structure.atoms[attach_id].pos)
                    attach_normal = structure.atoms[attach_to].normal
                    extracted_atoms = new_mof[attach_id:attach_id+1]
                    new_mof[attach_id:attach_id+1] = [None]
                    for _trial_rotation in range(rotations):
                        incoming_group = attachment.atoms_attached_to(attach_at, attach_towards, attach_normal)
                        for atom in incoming_group:
                            if not test_collision(atom, new_mof, structure.cell):
                                debug("Rotating group")
                                attach_normal = dot(rotation_about_angle(attach_towards, rotation_angle), attach_normal)
                                # Don't need to test more atoms
                                break
                        else:
                            # Fits, so add and move on
                            new_mof.extend(incoming_group)
                            break
                    else:
                        did_not_attach = this_group
                        break
                if did_not_attach:
                    break
            new_mof_name = "".join(new_mof_name)
            if did_not_attach:
                error("%i failed: %s from %s" % (idx, this_group, group_set))
                continue
            new_mof = [an_atom for an_atom in new_mof if an_atom is not None]
            info("%i: %s" % (idx, new_mof_name))
            job_name = structure.name
            with open('%s_func_%05i.xyz' % (job_name, idx), 'wb') as output_file:
                output_file.writelines(to_xyz(new_mof, name=new_mof_name))
            with open('%s_func_%05i.pdb' % (job_name, idx), 'wb') as output_file:
                output_file.writelines(to_pdb(new_mof, structure.cell, name=new_mof_name))


def random_replace(structure, groups, count=None, rotations=36):
    """
    Replace a random number of sites.

    """
    nsites = sum(len(x) for x in structure.attachments.values())
    if count is None:
        count = random.randint(1, nsites)
    elif count > nsites:
        warn("too many sites requested")
        count = nsites
    func_repr = [random.choice(groups.keys()) for _counter in range(count)] + [None]*(nsites - count)
    random.shuffle(func_repr)
    unique_name = hashlib.md5(str(func_repr)).hexdigest()
    new_mof_name = []
    new_mof = list(structure.atoms)
    for this_point, this_group in zip(chain(*[structure.attachments[x] for x in sorted(structure.attachments)]), func_repr):
        if this_group is None:
            new_mof_name.append("")
            continue
        else:
            new_mof_name.append(this_group)
        attachment = groups[this_group]
        attach_id = this_point[0]
        attach_to = this_point[1][0][1]
        attach_at = structure.atoms[attach_to].pos
        attach_towards = direction3d(attach_at, structure.atoms[attach_id].pos)
        attach_normal = structure.atoms[attach_to].normal
        extracted_atoms = new_mof[attach_id:attach_id+1]
        new_mof[attach_id:attach_id+1] = [None]
        for trial_rotation in range(rotations):
            incoming_group = attachment.atoms_attached_to(attach_at, attach_towards, attach_normal)
            for atom in incoming_group:
                if not test_collision(atom, new_mof, structure.cell):
                    debug("Rotating group")
                    attach_normal = dot(rotation_about_angle(attach_towards, random.random()*np.pi*2), attach_normal)
                    break
            else:
                new_mof.extend(incoming_group)
                break
        else:
            # this_point not valid
            did_not_attach = this_group
            error("Failed to generate: %s" % ".".join([x or "" for x in func_repr]))
            warn("Stopped after: %s" % ".".join(new_mof_name))

            return False

    new_mof = [an_atom for an_atom in new_mof if an_atom is not None]
    new_mof_name = "{" + ".".join(new_mof_name) + "}"
    info("Generated: %s" % new_mof_name)
    info("With unique name: %s" % unique_name)
    with open('random-%s.xyz' % unique_name, 'wb') as output_file:
        output_file.writelines(to_xyz(new_mof, name=new_mof_name))
    with open('random-%s.pdb' % unique_name, 'wb') as output_file:
        output_file.writelines(to_pdb(new_mof, structure.cell, name=new_mof_name))

    # completed sucessfully
    return True


def main():
    """
    Run the substitution for an input structure.

    """

    job_options = Options()
    job_name = job_options.get('job_name')

    input_structure = ModifiableStructure(job_name)
    input_structure.from_file(job_name,
                              job_options.get('initial_structure_format'),
                              job_options)
    #input_structure.from_cif("test_cifs/CALF21.cif")

    input_structure.gen_site_connection_table()
    input_structure.gen_normals()

    f_groups = FunctionalGroupLibrary()
    f_groups.from_file()

    all_combinations_replace(input_structure, f_groups)
#    for _make_some_randoms in range(30):
#        random_replace(input_structure, f_groups)


if __name__ == '__main__':

    main()
