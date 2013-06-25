#!/usr/bin/env python

"""
The basic Atom class that eveything is made up of.

"""

import re
from logging import warning, debug, error, info, critical

from numpy import identity, dot


class Atom(object):
    """Base atom object."""

    def __init__(self, at_type=False, pos=False, parent=None, **kwargs):
        """Accept arbritary kwargs as attributes."""
        if parent is not None:
            self._parent = parent
        self.type = at_type
        self.pos = pos
        self.charge = 0.0
        self.idx = False
        self.site = None
        self.mass = 0.0
        self.molecule = None
        self.uff_type = None
        self.is_fixed = False
        # Sets anything else specified as an attribute
        for key, val in kwargs.items():
            setattr(self, key, val)

    def __str__(self):
        return "%s %f %f %f" % tuple([self.type] + list(self.pos))

    def __repr__(self):
        return "Atom(%r,%r)" % (self.type, self.pos)

    def fpos(self, inv_cell):
        """Fractional position within a given cell."""
        return dot(inv_cell, self.pos).tolist()

    def ifpos(self, inv_cell):
        """In cell fractional position."""
        return [i % 1 for i in self.fpos(inv_cell)]

    def ipos(self, cell, inv_cell):
        """In cell cartesian position."""
        return dot(self.ifpos(inv_cell), cell)

    def from_cif(self, at_dict, cell, symmetry=None, idx=None):
        """Extract an atom description from dictionary of cif items."""
        self.site = at_dict['_atom_site_label']
        # type_symbol takes precedence but need not be specified
        self.type = at_dict.get('_atom_site_type_symbol', self.site)
        self.mass = WEIGHT[self.type]
        frac_pos = [ufloat(at_dict['_atom_site_fract_x']),
                    ufloat(at_dict['_atom_site_fract_y']),
                    ufloat(at_dict['_atom_site_fract_z'])]
        if symmetry is not None:
            frac_pos = symmetry.trans_frac(frac_pos)
        self.pos = dot(frac_pos, cell)
        if re.match('[0-9]', self.site) and idx is not None:
            debug("Site label may not be unique; appending index")
            self.site = "%s%i" % (self.site, idx)
        if '_atom_site_description' in at_dict:
            self.uff_type = at_dict['_atom_site_description']
        elif '_atom_type_description' in at_dict:
            self.uff_type = at_dict['_atom_type_description']
        if '_atom_type_partial_charge' in at_dict:
            self.charge = float(at_dict['_atom_type_partial_charge'])
        elif '_atom_type_parital_charge' in at_dict:
            #TODO(tdaff): remove for 2.0
            self.charge = float(at_dict['_atom_type_parital_charge'])

    def from_pdb(self, line, charges=False):
        """
        Parse the ATOM line from a pdb file.
        Occupancy field may be used to specify the charge as in a '.pqr' file.
        """
        # pdb is defined with fixed width fields rather than splitting
        self.idx = try_int(line[6:11])
        self.site = line[12:16].strip()
        self.molecule = try_int(line[22:26])
        at_pos = float(line[30:38]), float(line[38:46]), float(line[47:54])
        self.pos = at_pos
        self.type = line[76:78].strip()
        self.mass = WEIGHT[self.type]
        if charges:
            self.charge = float(line[54:60])

    def from_vasp(self, line, at_type=None, cell=identity(3)):
        """Set the atom data from vasp input. Only pass cell if fractional."""
        self.pos = dot([float(x) for x in line.split()[:3]], cell)
        if at_type is not None:
            self.type = at_type
            self.mass = WEIGHT[at_type]

    def from_siesta(self, line, cell):
        """Parse line from SIESTA.STRUCT_OUT file."""
        self.pos = dot([float(x) for x in line.split()[2:5]], cell)
        self.atomic_number = int(line.split()[1])
        self.mass = WEIGHT[self.type]

    def from_xyz(self, line):
        """Parse line from generic xyz file."""
        split_line = line.split()
        self.pos = [float(x) for x in split_line[1:4]]
        if len(split_line) > 4:
            try:
                self.charge = float(split_line[4])
            except ValueError:
                # Assume a comment so skip it
                pass
        self.type = line.split()[0]
        self.mass = WEIGHT[self.type]

    def translate(self, vec):
        """Move the atom by the given vector."""
        self.pos = [x + y for x, y in zip(self.pos, vec)]

    def get_atomic_number(self):
        """The atomic number for the element, or closest match."""
        name = self.type
        atomic_number = None
        while name:
            try:
                atomic_number = ATOMIC_NUMBER.index(name)
                break
            except ValueError:
                name = name[:-1]
        return atomic_number

    def set_atomic_number(self, value):
        """Set the atom type based on the atomic number."""
        self.type = ATOMIC_NUMBER[value]

    atomic_number = property(get_atomic_number, set_atomic_number)

    def get_fractional_coordinate(self):
        """Retrieve the fractional coordinates or calculate from the parent."""
        try:
            # Hopefully the attribute is just set
            return self._fractional
        except AttributeError:
            # Not set yet
            try:
                self._fractional = self.fpos(self._parent.cell.inverse)
                return self._fractional
            except AttributeError:
                return None

    def set_fractional_coordinate(self, value):
        """Set the position using the fractional coordinates."""
        fractional = [x % 1.0 for x in value]
        self._fractional = fractional
        self.pos = dot(fractional, self._parent.cell.cell)

    def del_fractional_coordinate(self):
        """Remove the fractional coordinate; run after updating cell."""
        del self._fractional

    fractional = property(get_fractional_coordinate, set_fractional_coordinate, del_fractional_coordinate)

    @property
    def element(self):
        """Guess the element from the type, fall back to type."""
        name = self.type
        while name:
            if name in ATOMIC_NUMBER:
                return name
            else:
                name = name[:-1]
        return self.type

    @property
    def vdw_radius(self):
        """Get the vdw radius from the UFF parameters."""
        return UFF[self.type][0]/2.0

    @property
    def covalent_radius(self):
        """Get the covalent radius from the library parameters."""
        if self.type == 'C' and self.uff_type:
            COVALENT_RADII[self.uff_type]
        return COVALENT_RADII[self.type]

    @property
    def is_metal(self):
        """Return True if element is in a predetermined set of metals."""
        return self.atomic_number in METALS
