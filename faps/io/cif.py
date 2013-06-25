"""

Reader for CIFs (Crystallographic information files).

Parses symmetry and some custom properties.

"""
import shlex
from logging import debug, info, warning, error

from faps.core.structure import Structure
from faps.core.symmetry import Symmetry
from faps.core.cell import Cell
from faps.util.text import strip_blanks, ufloat
from . import StructureReader

class CifReader(StructureReader):
    """
    CIF file reader with some symmetry and custom properties.

    """

    filetypes = ['*.cif', 'cif']

    @staticmethod
    def read_file(filename=None, string=None):
        """Genereate structure from a .cif file."""

        new_structure = Structure()

        if filename is not None:
            info("Reading positions from cif file: %s" % filename)
            with open(filename) as temp_file:
                cif_file = temp_file.readlines()
        elif string is not None:
            info("Positions from cif string")
            cif_file = string.splitlines()
        else:
            error("No source for cif file")
        cif_file = strip_blanks(cif_file)
        cell = Cell()
        params = [None, None, None, None, None, None]
        atoms = []
        cif_bonds = {}
        symmetry = Symmetry()
        loops = []
        idx = 0
        while idx < len(cif_file):
            # line text needs to be manageable; cif guidelines can be
            # permissive
            # Can no longer just check for _underscores in lines
            # as UFF types can have them and mess up parsing
            line = cif_file[idx].lower().strip()
            if '_cell_length_a' in line:
                params[0] = ufloat(line.split()[1])
            elif '_cell_length_b' in line:
                params[1] = ufloat(line.split()[1])
            elif '_cell_length_c' in line:
                params[2] = ufloat(line.split()[1])
            elif '_cell_angle_alpha' in line:
                params[3] = ufloat(line.split()[1])
            elif '_cell_angle_beta' in line:
                params[4] = ufloat(line.split()[1])
            elif '_cell_angle_gamma' in line:
                params[5] = ufloat(line.split()[1])
            elif '_symmetry_space_group_name_h-m' in line:
                Symmetry.space_group = line.split()[1]
            elif 'loop_' in line:
                # loops for _atom_site, _symmetry and _geom
                heads = []
                body = []
                while line.startswith('_') or 'loop_' in line:
                    # must keep the loop_ line as this can still contain headers
                    heads.extend(line.split())
                    idx += 1
                    # don't lower these to keep atomic symbols
                    line = cif_file[idx].strip()
                while idx < len(cif_file) and not line.startswith('_') and not 'loop_' in line:
                    # shlex keeps 'quoted items' as one
                    # Some cifs seem to have primed atom symbols
                    # posix=False should help
                    # using .shlex instead of .split works with '#' comments too
                    split_line = shlex.shlex(line, posix=False)
                    split_line.whitespace_split = True
                    split_line = list(split_line)
                    body.extend([x.strip("'").strip('"') for x in split_line])
                    idx += 1
                    try:
                        line = cif_file[idx]
                    except IndexError:
                        line = ''
                if 'loop_' in heads:
                    heads.remove('loop_')
                loops.append((heads, body))
                continue
            idx += 1

        # cell first
        if all(params):
            cell.params = params
            new_structure.cell = cell
        else:
            error("No cell or incomplete cell found in cif file")

        # parse loop contents
        for heads, body in loops:
            if '_atom_site_fract_x' in heads:
                while body:
                    atoms.append(dict(zip(heads, body)))
                    body = body[len(heads):]
            if '_symmetry_equiv_pos_as_xyz' in heads:
                while body:
                    sym_dict = dict(zip(heads, body))
                    symmetry.append(
                        Symmetry(sym_dict['_symmetry_equiv_pos_as_xyz']))
                    body = body[len(heads):]
            if '_ccdc_geom_bond_type' in heads:
                while body:
                    bond_dict = dict(zip(heads, body))
                    # bond is sorted so there are no duplicates
                    # and tuple so it can be hashed
                    bond = (bond_dict['_geom_bond_atom_site_label_1'],
                            bond_dict['_geom_bond_atom_site_label_2'])
                    bond = tuple(sorted(bond))
                    # bond distance and type defualts to None if not specified
                    distance = bond_dict.get('_geom_bond_distance')
                    if distance is not None:
                        distance = float(distance)
                    bond_type = bond_dict.get('_ccdc_geom_bond_type')
                    cif_bonds[bond] = (distance, bond_type)
                    body = body[len(heads):]

        if not symmetry:
            debug('No symmetry found; assuming identity only')
            symmetry = [Symmetry('x,y,z')]

        duplicate_tolerance = 0.2  # Angstroms
        symmetry.expand(atoms, duplicate_tolerance)

        new_structure.atoms = atoms

        bonds = {}
        # TODO(tdaff): this works for the one tested MOF; 0.1 was not enough
        # only check for bonds that are too long, not too short.
        # FIXME(tdaff): bonds with two distances but same IDX
        bond_tolerence = 0.25
        # Assign bonds by index
        for bond, bond_data in cif_bonds.items():
            for first_index, first_atom in enumerate(atoms):
                if first_atom.site == bond[0]:
                    for second_index, second_atom in enumerate(atoms):
                        if second_atom is first_atom:
                            continue
                        elif second_atom.site == bond[1]:
                            # TODO(tdaff): symmetry implementation for cif bonding
                            distance = min_distance(first_atom, second_atom)
                            bond_dist = bond_data[0]
                            if bond_dist is None:
                                bond_dist = first_atom.covalent_radius + second_atom.covalent_radius
                            if distance < (bond_dist + bond_tolerence):
                                # use the sorted index as bonds between the
                                # same type are doubly specified
                                bond_id = tuple(sorted((first_index, second_index)))
                                bonds[bond_id] = CCDC_BOND_ORDERS[bond_data[1]]
                                if first_atom.is_metal or second_atom.is_metal:
                                    first_atom.is_fixed = True
                                    second_atom.is_fixed = True

        new_structure.bonds = bonds
        new_structure.symmetry = symmetry

        return new_structure