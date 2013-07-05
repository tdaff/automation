"""
The core Structure is the central basis for storing atomic configurations
and everything attached to them, including the route to derive the
coordinates and the attributes.

"""


from copy import copy
from itertools import count
from logging import warning, debug, info, critical
from os import path

from numpy import dot


class Structure(object):
    """
    The current state of the structure.

    Structure holds data on atomic postitions, it should also optionally 
    define unique sets of attributes for each configuration.
    Internal energy units are kcal/mol.

    """

    def __init__(self, drepr=None):
        """Instance an empty container with nothing in it."""
        super(Structure, self).__setattr__('_attribute_register',set())
        self.name = None
        self.cell = None
        self.atoms = []
        self.properties = {}

    def __setattr__(self, name, value):
        self._attribute_register.add(name)
        return super(Structure, self).__setattr__(name, value)

    def __getattribute__(self, name):
        print("getarrt")
        return super(Structure, self).__getattribute__(name)

    def remove_duplicates(self, tolerance=0.02):
        """Find overlapping atoms and remove them."""
        unique_atoms = []
        found_atoms = []
        for atom in self.atoms:
            for unique_atom in unique_atoms:
                if atom.type != unique_atom.type:
                    continue
                elif min_distance(atom, unique_atom) < tolerance:
                    break
            # else excutes when not found here
            else:
                unique_atoms.append(atom)
        debug("Found %i unique atoms in %i" % (len(unique_atoms), self.natoms))
        self.atoms = unique_atoms

    def check_close_contacts(self, absolute=1.0, covalent=None):
        """
        Check for atoms that are too close. Specify either an absolute distance
        in Angstrom or a scale factor for the sum of covalent radii. If a
        covalent factor is specified it will take priority over an absolute
        distance. Return True if close contacts found, else return False.
        """
        close_contact_found = False
        for atom_idx, atom in enumerate(self.atoms):
            for other_idx, other in enumerate(self.atoms):
                if other_idx >= atom_idx:
                    # short circuit half the calculations
                    # Can we do combinations with idx in 2.7
                    break
                if covalent is not None:
                    tolerance = covalent * (atom.covalent_radius +
                                            other.covalent_radius)
                else:
                    tolerance = absolute
                if min_distance(atom, other) < tolerance:
                    bond_ids = tuple(sorted([atom_idx, other_idx]))
                    if bond_ids not in self.bonds:
                        warning("Close atoms: %s(%i) and %s(%i)" %
                                (atom.site, atom_idx, other.site, other_idx))
                        close_contact_found = True

        return close_contact_found

    def bond_length_check(self, too_long=1.25, too_short=0.7):
        """
        Check if all bonds fall within a sensible range of scale factors
        of the sum of the covalent radii. Return True if bad bonds are found,
        otherwise False.

        """
        bad_bonds = False
        for bond in self.bonds:
            atom = self.atoms[bond[0]]
            other = self.atoms[bond[1]]
            distance = min_distance(atom, other)
            bond_dist = (atom.covalent_radius + other.covalent_radius)
            if distance > bond_dist * too_long:
                warning("Long bond found: %s(%i) and %s(%i) = %.2f A" %
                        (atom.site, bond[0], other.site, bond[1], distance))
                bad_bonds = True
            elif distance < bond_dist * too_short:
                warning("Short bond found: %s(%i) and %s(%i) = %.2f A" %
                        (atom.site, bond[0], other.site, bond[1], distance))
                bad_bonds = True

        return bad_bonds

    def gen_supercell(self, options):
        """Cacluate the smallest satisfactory supercell and set attribute."""
        config_supercell = options.gettuple('mc_supercell', int)
        config_cutoff = options.getfloat('mc_cutoff')
        if config_cutoff < 12:
            warning("Simulation is using a very small cutoff! I hope you "
                 "know, what you are doing!")
        minimum_supercell = self.cell.minimum_supercell(config_cutoff)
        supercell = tuple(max(i, j)
                          for i, j in zip(config_supercell, minimum_supercell))
        self.gcmc_supercell = supercell
        info("%s supercell requested in config" % str(config_supercell))
        info("%s minimum supercell for a %.1f cutoff" %
             (str(minimum_supercell), config_cutoff))
        info("Constructing %s supercell for gcmc." % str(supercell))

    def supercell(self, scale):
        """
        Iterate over all the atoms of supercell where scale is an integer
        to scale uniformly or triplet with scale factors for each direction.

        """
        # Beware supercells larger than 2147483647 are not supported in
        # python 2
        if isinstance(scale, int):
            scale = (scale, scale, scale)
        for x_super in range(scale[0]):
            for y_super in range(scale[1]):
                for z_super in range(scale[2]):
                    offset = dot((x_super, y_super, z_super), self.cell.cell)
                    for atom in self.atoms:
                        newatom = copy(atom)
                        newatom.translate(offset)
                        yield newatom

    def order_by_types(self):
        """Sort the atoms alphabetically and group them."""
        self.atoms.sort(key=lambda x: (x.type, x.site))

    def gen_neighbour_list(self, force=False):
        """All atom pair distances."""
        # This can be expensive so skip if already calcualted
        if not force:
            for atom in self.atoms:
                if not hasattr(atom, 'neighbours'):
                    break
            else:
                # finished loop over all atoms
                debug("Neighbour list already calculated")
                return

        debug("Calculating neighbour list.")

        cell = self.cell.cell
        inv_cell = self.cell.inverse
        cpositions = [atom.ipos(cell, inv_cell) for atom in self.atoms]
        fpositions = [atom.ifpos(inv_cell) for atom in self.atoms]
        cell = cell.tolist()

        # loop over all pairs to find minimum periodic distances
        for atom, a_cpos, a_fpos in zip(self.atoms, cpositions, fpositions):
            neighbours = []
            for o_idx, o_cpos, o_fpos in zip(count(), cpositions, fpositions):
                sep = min_dist(a_cpos, a_fpos, o_cpos, o_fpos, cell)
                neighbours.append((sep, o_idx))
            # First one is self == 0
            # save in incresaing distance order
            atom.neighbours = sorted(neighbours)[1:]

        ## Neighbourlist printed in VASP style
        #for idx, atom in enumerate(self.atoms):
        #    print("%4i" % (idx+1) +
        #          "%7.3f%7.3f%7.3f" % tuple(atom.ifpos(inv_cell)) +
        #          "-" +
        #          "".join("%4i%5.2f" % (y+1, x) for x, y in atom.neighbours if x<2.5))

    def surface_area(self, probe=None, value=None, delete=False):
        """
        Helper:
          Return all {probe:area} if no arguments given
          Return the area or None for a given probe
          Set area if value given
          Delete value if delete is True
        Areas in A^2
        """
        surface_areas = self.properties.get('surface_area', {})
        if value is not None:
            surface_areas[probe] = value
            self.properties['surface_area'] = surface_areas
        elif delete:
            # Set it to None to avoid KeyErrors
            surface_areas[probe] = None
            del surface_areas[probe]
            self.properties['surface_area'] = surface_areas
        elif probe is not None:
            return surface_areas.get(probe, None)
        else:
            return surface_areas

    def sub_property(self, name, probe=None, value=None, delete=False):
        """
        Helper:
          Return all {probe:value} if no arguments given
          Return the value or None for a given probe
          Set area if value given
          Delete value if delete is True
        Units are based on Angstrom
        """
        property_data = self.properties.get(name, {})
        if value is not None:
            property_data[probe] = value
            self.properties[name] = property_data
        elif delete:
            # Set it to None to avoid KeyErrors
            property_data[probe] = None
            del property_data[probe]
            self.properties[name] = property_data
        elif probe is not None:
            return property_data.get(probe, None)
        else:
            return property_data

    @property
    def types(self):
        """Ordered list of atom types."""
        return [atom.type for atom in self.atoms]

    @property
    def atomic_numbers(self):
        """Ordered list of atomic numbers."""
        return [atom.atomic_number for atom in self.atoms]

    @property
    def weight(self):
        """Unit cell weight."""
        return sum([atom.mass for atom in self.atoms])

    @property
    def volume(self):
        """Unit cell volume."""
        return self.cell.volume

    @property
    def natoms(self):
        """Number of atoms in the unit cell."""
        return len(self.atoms)

    @property
    def symmetry_tree(self):
        """Tree of atoms that are symmetrically equivalent."""
        tree = {}
        for atom_id, atom in enumerate(self.atoms):
            if atom.site in tree:
                tree[atom.site].append(atom_id)
            else:
                tree[atom.site] = [atom_id]
        if len(tree) == 1 and None in tree:
            return dict((i, [i]) for i in range(self.natoms))
        else:
            return tree

    def get_gcmc_supercell(self):
        """Supercell used for gcmc."""
        return self.properties.get('supercell', (1, 1, 1))

    def set_gcmc_supercell(self, value):
        """Set the supercell property for the structure."""
        self.properties['supercell'] = value

    gcmc_supercell = property(get_gcmc_supercell, set_gcmc_supercell)

    @property
    def to_dict(self):
        """
        Dictionary representation of structure for document interchange.
        """
        drepr = {"@module": self.__class__.__module__,
                 "@class": self.__class__.__name__,
                 "cell": self.cell.to_dict or None}
        return drepr

    @staticmethod
    def from_dict(drepr):
        return Structure(drepr=drepr)


if __name__ == "__main__":
    lol = Structure()
    lol.cow = 12
    print(lol)
    print(lol.__dict__)
    print(lol.to_dict)
