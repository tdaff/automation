#!/usr/bin/env python

"""
The core Structure as the central basis for storing atomic configurations.

"""

import code
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import glob
import os
import pickle
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import textwrap
import time
from copy import copy
from itertools import count
from logging import warning, debug, error, info, critical
from math import ceil, log
from os import path

import numpy as np
from numpy import pi, cos, sin, sqrt, arccos
from numpy import array, identity, dot, cross
from numpy.linalg import norm

from config import Options
from elements import WEIGHT, ATOMIC_NUMBER, UFF, VASP_PSEUDO_PREF
from elements import CCDC_BOND_ORDERS, GULP_BOND_ORDERS, METALS
from elements import COVALENT_RADII, UFF_FULL, QEQ_PARAMS
from job_handler import JobHandler
from logo import LOGO

# Global constants
DEG2RAD = pi / 180.0
BOHR2ANG = 0.52917720859
EV2KCAL = 23.060542301389
NAVOGADRO = 6.02214129E23
INFINITY = float('inf')
KCAL_TO_KJ = 4.1868  # Steam tables from openbabel

FASTMC_DEFAULT_GRID_SPACING = 0.1


class Structure(object):
    """
    The current state of the structure; update as the calculations proceed.

    Structure provides methods to produce input files for and take output from
    various computational chemistry packages but needs to be told what to do.
    Internal energy units are kcal/mol.

    Methods are grouped:
    * Initial structure parsers
    * Output file parsers to update structure
    * Input file generation
    * Internal manipulation methods

    """
    # TODO: dft energy?
    def __init__(self, name):
        """Just instance an empty structure initially."""
        self.name = name
        self.cell = Cell()
        self.atoms = []
        self.esp = None
        self.dft_energy = 0.0
        self.guests = []
        self.properties = {}
        self.space_group = None

    def from_file(self, basename, filetype, defaults):
        """Select the correct file parser."""
        if filetype in ['sql', 'sqlite']:
            # Makeshift selection method to select a mof
            # jobname is dbname.type.structure_id
            from backend.sql import AlchemyBackend
            # [db_name, sym or free, identity]
            db_params = self.name.split('.')
            reader = AlchemyBackend(db_params[0])
            cif_string = reader.start_cif(db_params[1], int(db_params[2]))
            self.from_cif(string=cif_string)
        elif filetype.lower() in ['pdb']:
            self.from_pdb(basename + '.' + filetype)
        elif filetype.lower() in ['pqr']:
            # Look for a pqr or just a pdb wih charges
            listdir = os.listdir('.')
            if (basename + '.pqr') in listdir:
                self.from_pdb(basename + '.pqr', charges=True)
            else:
                self.from_pdb(basename + '.pdb', charges=True)
        elif filetype.lower() in ['vasp', 'poscar', 'contcar']:
            listdir = os.listdir('.')
            test_files = [
                basename + '.contcar', basename + '.CONTCAR', 'CONTCAR',
                basename + '.poscar', basename + '.POSCAR', 'POSCAR']
            for filename in test_files:
                if filename in listdir:
                    self.from_vasp(filename)
                    break
        elif filetype.lower() in ['cif']:
            self.from_cif(basename + '.' + filetype)
        elif filetype.lower() in ['xyz']:
            cell = defaults.gettuple('default_cell', float)
            self.from_xyz(basename + '.' + filetype, cell=cell)
        else:
            error("Unknown filetype %s" % filetype)

    def update_pos(self, opt_code):
        """Select the method for updating atomic positions."""
        opt_path = path.join('faps_%s_%s' % (self.name, opt_code))
        info("Updating positions from %s" % opt_code)
        if opt_code == 'vasp':
            self.from_vasp(path.join(opt_path, 'CONTCAR'), update=True)
        elif opt_code == 'siesta':
            self.from_siesta(path.join(opt_path, '%s.STRUCT_OUT' % self.name))
        elif opt_code == 'gulp':
            opt_path = "%s_opt" % opt_path
            self.optimisation_output = validate_gulp_output(
                path.join(opt_path, 'faps-%s.out' % self.name))
            self.from_gulp_output(path.join(opt_path, '%s.grs' % self.name))
        else:
            error("Unknown positions to import %s" % opt_code)

    def update_charges(self, charge_method, options=None):
        """Select the method for updating charges."""
        charge_path = path.join('faps_%s_%s' % (self.name, charge_method))
        if charge_method == 'repeat':
            info("Updating charges from repeat")
            self.charges_from_repeat(
                path.join(charge_path, 'faps-%s.out' % self.name),
                options.getbool('symmetry'))
            # Cleanup of REPEAT files
            unneeded_files = options.gettuple('repeat_delete_files')
            remove_files(unneeded_files, charge_path)
            keep_files = options.gettuple('repeat_compress_files')
            compress_files(keep_files, charge_path)
        elif charge_method == 'gulp':
            info("Updating charges from GULP QEq")
            self.charges_from_gulp(
                path.join(charge_path, 'faps-%s.out' % self.name))
        elif charge_method == 'egulp':
            info("Updating charges from EGULP QEq")
            self.charges_from_egulp(path.join(charge_path, 'charges.dat'))
        else:
            error("Unknown charge method to import %s" % charge_method)

    def update_gcmc(self, tp_point, options):
        """Select the source for GCMC results and import."""
        gcmc_code = options.get('mc_code')
        gcmc_path = path.join('faps_%s_%s' % (self.name, gcmc_code))
        # Runs in subdirectories
        tp_path = path.join(gcmc_path, 'T%s' % tp_point[0] +
                               ''.join(['P%.2f' % x for x in tp_point[1]]))
        if gcmc_code == 'fastmc':
            info("Importing results from FastGCMC")
            self.fastmc_postproc(tp_path, tp_point, options)
        else:
            error("Unknown gcmc method to import %s" % gcmc_code)

    def from_pdb(self, filename, charges=False):
        """Read an initial structure from a pdb file."""
        info("Reading positions from pdb file: %s" % filename)
        filetemp = open(filename)
        pdb_file = filetemp.readlines()
        filetemp.close()
        # Build a local list before setting attribute
        newatoms = []
        for line in pdb_file:
            lline = line.lower()
            if lline.startswith('cryst1'):
                self.cell.from_pdb(line)
                self.space_group = line[55:56]
            elif lline.startswith('atom') or lline.startswith('hetatm'):
                newatom = Atom()
                newatom.from_pdb(line, charges=charges)
                newatoms.append(newatom)

        self.atoms = newatoms
        self.order_by_types()

    def from_cif(self, filename=None, string=None):
        """Genereate structure from a .cif file."""
        if filename is not None:
            info("Reading positions from cif file: %s" % filename)
            filetemp = open(filename)
            cif_file = filetemp.readlines()
            filetemp.close()
        elif string is not None:
            info("Positions from cif string")
            cif_file = string.splitlines()
        else:
            error("No source for cif file")
        cif_file = strip_blanks(cif_file)
        params = [None, None, None, None, None, None]
        atoms = []
        cif_bonds = {}
        symmetry = []
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
                self.space_group = line.split()[1]
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
        if np.all(params):
            self.cell.params = params
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

        newatoms = []
        for site_idx, atom in enumerate(atoms):
            for sym_op in symmetry:
                newatom = Atom(parent=self)
                newatom.from_cif(atom, self.cell.cell, sym_op, site_idx)
                newatoms.append(newatom)

        self.atoms = newatoms

        if len(symmetry) > 1:
            # can skip if just identity operation as it's slow for big systems
            # Some of pete's symmetrised mofs need a higher tolerence
            duplicate_tolerance = 0.2  # Angstroms
            self.remove_duplicates(duplicate_tolerance)
        self.order_by_types()

        bonds = {}
        # TODO(tdaff): this works for the one tested MOF; 0.1 was not enough
        # only check for bonds that are too long, not too short.
        bond_tolerence = 0.25
        # Assign bonds by index
        for bond, bond_data in cif_bonds.items():
            for first_index, first_atom in enumerate(self.atoms):
                if first_atom.site == bond[0]:
                    for second_index, second_atom in enumerate(self.atoms):
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

        self.bonds = bonds
        self.symmetry = symmetry

    def from_vasp(self, filename='CONTCAR', update=False):
        """Read a structure from a vasp [POS,CONT]CAR file."""
        info("Reading positions from vasp file: %s" % filename)
        filetemp = open(filename)
        contcar = filetemp.readlines()
        filetemp.close()
        atom_list = []
        atom_types = []
        scale = float(contcar[1])
        self.cell.from_lines(contcar[2:5], scale)
        if contcar[5].split()[0].isalpha():
            # vasp 5 with atom names
            atom_types = contcar[5].split()
            del contcar[5]
        poscar_counts = [int(x) for x in contcar[5].split()]
        natoms = sum(poscar_counts)
        if contcar[6].strip()[0].lower() in "s":
            # 's'elective dynamics line; we don't care
            del contcar[6]

        # mcell converts frac -> cart if necessary and scales
        if contcar[6].strip()[0].lower() in "ck":
            mcell = identity(3) * scale
        else:
            mcell = self.cell.cell

        # parsing positions
        if update:
            for atom, at_line in zip(self.atoms, contcar[7:7+natoms]):
                atom.from_vasp(at_line, cell=mcell)
        elif not atom_types:
            critical("Will not extract structure from older vasp files")
        else:
            line_idx = 6
            for at_type, at_count in zip(atom_types, poscar_counts):
                for _atom_idx in range(at_count):
                    line_idx += 1
                    this_atom = Atom()
                    this_atom.from_vasp(contcar[line_idx], at_type, mcell)
                    atom_list.append(this_atom)
            self.atoms = atom_list
            self.order_by_types()

    def from_siesta(self, filename):
        """Update the structure from the siesta output."""
        info("Updating positions from file: %s" % filename)
        filetemp = open(filename)
        struct_out = filetemp.readlines()
        filetemp.close()
        self.cell.from_lines(struct_out[:3])
        for atom, line in zip(self.atoms, struct_out[4:]):
            atom.from_siesta(line, self.cell.cell)

    def from_gulp_output(self, filename):
        """Update the structure from the gulp optimisation output."""
        info("Updating positions from file: %s" % filename)
        grs_out = open(filename)
        for line in grs_out:
            if line.strip() == 'cell':
                params = tuple(float(x) for x in grs_out.next().split()[:6])
                self.cell.params = params
            elif line.strip() == 'fractional':
                cell = self.cell.cell
                for atom, atom_line in zip(self.atoms, grs_out):
                    atom.pos = dot([gfloat(x) for x in atom_line.split()[2:5]], cell)
                    # FIXME(tdaff) fractionals need to change automatically
                    # when the position gets updated (or the cell!)
                    del atom.fractional

        # Make sure everything is good from here
        if self.check_close_contacts(covalent=1.0):
            warning("Structure might have atom overlap, check gulp output!")
            self.bad_structure = True

        if self.bond_length_check():
            warning("Structure might have strained bonds, check gulp output!")
            self.bad_structure = True


    def from_xyz(self, filename, update=False, cell=None):
        """Read a structure from an file."""
        info("Reading positions from xyz file: %s" % filename)
        filetemp = open(filename)
        xyz_file = filetemp.readlines()
        filetemp.close()
        # Setting the cell
        if len(cell) == 6:
            self.cell.params = cell
        elif len(cell) == 9:
            self.cell.cell = array(cell).reshape((3, 3))
        elif cell is not None:
            error("Invalid cell specification %s" % str(cell))
        # Build a local list before setting attribute
        newatoms = []
        natoms = int(xyz_file[0])
        self.properties['header'] = xyz_file[1].strip()
        for line in xyz_file[2:2+natoms]:
            newatom = Atom()
            newatom.from_xyz(line)
            newatoms.append(newatom)
        if update:
            if natoms != self.natoms:
                critical("Incorrect number of atoms to update")
                terminate(96)
            for atom, newatom in zip(self.atoms, newatoms):
                if atom.type != newatom.type:
                    error("Atom order may have changed")
                atom.pos = newatom.pos
        else:
            self.atoms = newatoms
            self.order_by_types()

    def charges_from_repeat(self, filename, symmetry=False):
        """Parse charges and update structure."""
        info("Getting charges from file: %s" % filename)
        charges = []
        filetemp = open(filename)
        for line in filetemp:
            if line.startswith(" Charge"):
                line = line.split()
                # index, type, charge
                charges.append((int(line[1]), int(line[4]), float(line[6])))
            if "Error" in line:
                if float(line.split()[-1]) > 0.6:
                    warning("Error in repeat charges is very high - check cube!")
        filetemp.close()
        if symmetry:
            tree = self.symmetry_tree
            if len(charges) != len(tree):
                error("Incorrect number of charge sets; check REPEAT output")
                terminate(97)
            for symm, charge in zip(sorted(tree.items()), charges):
                for at_idx in symm[1]:
                    self.atoms[at_idx].charge = charge[2]
        else:
            if len(charges) != len(self.atoms):
                error("Incorrect number of charges; check REPEAT output")
                terminate(90)
            for atom, charge in zip(self.atoms, charges):
                atom.charge = charge[2]

    def charges_from_gulp(self, filename):
        """Parse QEq charges from GULP output."""
        info("Getting charges from file: %s" % filename)
        filetemp = open(filename)
        gout = filetemp.readlines()
        filetemp.close()
        start_line = None
        failures = 0
        for line_num, line in enumerate(gout):
            if '  Final charges from QEq :' in line:
                # Take the last set of charges
                start_line = line_num + 7
            elif 'Failed to converge' in line:
                # This will occur twice in the output for complete failure
                failures += 1
        if start_line is None:
            error("Final charges not found in gulp output")
            terminate(184)
        elif failures > 1:
            warning("Gulp charges may not be converged")
        for atom, chg_line in zip(self.atoms, gout[start_line:]):
            atom.charge = float(chg_line.split()[2])

    def charges_from_egulp(self, filename):
        """Parse QEq charges from EGULP output."""
        info("Getting charges from file: %s" % filename)
        filetemp = open(filename)
        gout = filetemp.readlines()
        filetemp.close()
        # charges.dat file only has list of charges in it
        for atom, chg_line in zip(self.atoms, gout):
            atom.charge = float(chg_line.split()[2])
            if atom.charge != atom.charge:
                # These are 'nan'; should not continue or fastmc will mess up
                error("Egulp charges did not converge, check structure")
                terminate(107)
            elif abs(atom.charge) == INFINITY:
                error("Egulp gave infinite charges, check structure")
                terminate(108)
            elif abs(atom.charge) > 10:
                warning("Very high charge from egulp: %s %f"
                        (atom.site, atom.charge))


    def to_vasp(self, options):
        """Return a vasp5 poscar as a list of lines."""
        optim_h = options.getbool('optim_h')
        optim_all = options.getbool('optim_all')
        types = [atom.type for atom in self.atoms]
        ordered_types = unique(types)
        poscar = ["%s\n" % self.name[:80],
                  " 1.0\n"]
        # Vasp does 16 dp but we get rounding errors (eg cubic) use 14
        poscar.extend(self.cell.to_vector_strings(fmt="%23.14f"))
        poscar.append("".join("%5s" % x for x in ordered_types) + "\n")
        poscar.append("".join("%6i" % types.count(x)
                              for x in ordered_types) + "\n")
        # We always have the T or F so turn on selective dynamics for
        # fixed pos variable cell
        poscar.extend(["Selective dynamics\n", "Cartesian\n"])
        if optim_all:
            info("Optimizing all atom positions")
            fix_h = 'T'
            fix_all = 'T'
        elif optim_h:
            info("Optimizing hydrogen positions")
            fix_h = 'T'
            fix_all = 'F'
        else:
            info("All atom positions fixed")
            fix_h = 'F'
            fix_all = 'F'

        # assume atoms are ordered
        # Atom lines differ from vasp to ensure spaces between numbers
        # with <-10 values
        for atom in self.atoms:
            if atom.type == "H":
                poscar.append("%20.15f %19.15f %19.15f" % tuple(atom.pos) +
                              "%4s%4s%4s\n" % (fix_h, fix_h, fix_h))
            else:
                poscar.append("%20.15f %19.15f %19.15f" % tuple(atom.pos) +
                              "%4s%4s%4s\n" % (fix_all, fix_all, fix_all))
        return poscar

    def to_siesta(self, options):
        """Return a siesta input file as a list of lines."""
        job_name = options.get('job_name')
        siesta_accuracy = options.get("siesta_accuracy").lower()
        if siesta_accuracy in ['high']:
            info("Using 'high' accuracy siesta settings")
            basis = ('DZP', 100, 200)
        elif siesta_accuracy in ['low']:
            info("Using 'low' accuracy siesta settings")
            basis = ('SZ', 200, 100)
        else:
            info("Using default siesta accuracy settings")
            basis = ('DZ', 150, 150)

        u_atoms = unique(self.atoms, key=lambda x: x.type)
        u_types = unique(self.types)
        fdf = ([
            "SystemName %s\n" % job_name,
            "SystemLabel %s\n" % job_name,
            "\n",
            "NumberOfAtoms %i\n" % len(self.atoms),
            "NumberOfSpecies %i\n" % len(u_atoms),
            "\n",
            "SaveElectrostaticPotential .true.\n",
            "WriteMullikenPop 1\n"
            "WriteXML F\n",
            "\n",
            "# Accuracy bits\n",
            "\n",
            "PAO.BasisSize %s\n" % basis[0],
            "PAO.EnergyShift %i meV\n" % basis[1],
            "MeshCutoff %i Ry\n" % basis[2],
            "\n",
            "MaxSCFIterations 100\n",
            "XC.Functional GGA\n",
            "XC.Authors PBE\n",
            "SolutionMethod diagon\n",
            "ElectronicTemperature 25 K\n",
            "DM.NumberPulay 5\n",
            "DM.MixingWeight 0.05\n",
            "\n",
            "%block ChemicalSpeciesLabel\n"] + [
            "%6i %6i %6s\n" % ((idx + 1), atom.atomic_number, atom.type)
                for idx, atom in enumerate(u_atoms)] + [
            "%endblock ChemicalSpeciesLabel\n",
            "\n",
            "LatticeConstant 1 Ang\n"
            "%block LatticeVectors\n"] +
            self.cell.to_vector_strings() + [
            "%endblock LatticeVectors\n",
            "\n",
            "AtomicCoordinatesFormat   Ang\n",
            "\n",
            "%block AtomicCoordinatesAndAtomicSpecies\n"] + [
            "%12.8f %12.8f %12.8f" % tuple(atom.pos) + "%6i\n" %
                (u_types.index(atom.type) + 1) for atom in self.atoms] + [
            "%endblock AtomicCoordinatesAndAtomicSpecies\n"])
        if "Br" in u_types:
            fdf.extend(["\n%block PS.lmax\n",
                        "   Br 2\n",
                        "%endblock PS.lmax\n"])

        if "In" in u_types:
            # semicore states need to be accounted for
            # TODO(tdaff): selective polarisation on the basis functions
            fdf.extend(["\n%block PAO.Basis\n",
                        "In   3\n",
                        "  n=5  0  2\n",
                        "    0.0 0.0\n",
                        "  n=5  1  2  P\n",
                        "    0.0 0.0\n",
                        "  n=4  2  2\n",
                        "    0.0 0.0\n",
                        "%endblock PAO.Basis\n"])

        optim_h = options.getbool('optim_h')
        optim_all = options.getbool('optim_all')
        optim_cell = options.getbool('optim_cell')

        constraint = []

        if optim_all:
            info("Optimizing all atom positions")
            fdf.append("\nMD.TypeOfRun  CG\n")
            fdf.append("MD.NumCGSteps  %i\n" % 800)
        elif optim_h and "H" in self.types:
            info("Optimizing hydrogen positions")
            fdf.append("\nMD.TypeOfRun  CG\n")
            fdf.append("MD.NumCGSteps  %i\n" % 300)
            constraint = ["%i" % (idx+1)
                          for idx, species in enumerate(self.types)
                          if species not in ["H"]]
        elif optim_cell:
            fdf.append("\nMD.TypeOfRun  CG\n")
            fdf.append("MD.NumCGSteps  %i\n" % 300)
            constraint = ["%i" % (idx+1) for idx in range(self.natoms)]

        if optim_cell:
            info("Cell vectors will be optimized")
            fdf.append("MD.VariableCell .true.\n")
        if constraint:
            constraint = textwrap.fill(" ".join(constraint),
                                       initial_indent='position ',
                                       subsequent_indent='position ')
            fdf.extend(["\n%block GeometryConstraints\n",
                        constraint, "\n",
                        "%endblock GeometryConstraints\n"])

        return fdf

    def to_gulp(self, qeq_fit=False, optimise=False, terse=False, qeq_dict={}):
        """Return a GULP file to use for the QEq charges."""
        if qeq_fit:
            keywords = "fitting bulk_noopt qeq\n"
        elif optimise:
            return self.to_gulp_optimise(terse=terse)
        else:
            keywords = "single conp qeq\n"

        gin_file = [
            "# \n# Keywords:\n# \n",
            keywords,
            "# \n# Options:\n# \n",
            "# Initial file written by faps\n"
            "name %s\n" % self.name,
            "vectors\n"] + self.cell.to_vector_strings() + [
            "cartesian\n"]
        if qeq_fit:
            for atom in self.atoms:
                gin_file.extend(["%-5s core " % atom.type,
                                 "%14.7f %14.7f %14.7f " % tuple(atom.pos),
                                 "%14.7f\n" % atom.charge])
            gin_file.append("\n# \n# Fitting variables:\n# \nobservables\n")
            for idx, atom in enumerate(self.atoms):
                gin_file.append("monopoleq\n%5s%11.6f\n" % (idx+1, atom.charge))
            gin_file.append("end\nqelectronegativity\n")
            for at_type in unique(self.types):
                gin_file.extend(
                    ["%-5s" % at_type,
                     "%9.5f %9.5f %9.5f 1 1 0\n" % QEQ_PARAMS[at_type]])
        else:
            for atom in self.atoms:
                gin_file.extend(["%-5s core " % atom.type,
                                 "%14.7f %14.7f %14.7f\n" % tuple(atom.pos)])
            if qeq_dict:
                unique_types = unique(self.types)
                gin_file.append('\nqelectronegativity\n')
                for atom_type, params in qeq_dict.items():
                    if atom_type in unique_types:
                        gin_file.append('%-4s %f %f\n' % (atom_type, params[0], params[1]))

        gin_file.append("\ndump every %s.grs\nprint 1\n" % self.name)
        return gin_file


    def to_gulp_optimise(self, terse=False):
        """Return a GULP file to optimise with UFF."""
        # all atoms of the same forcefield type must have the same label
        # gulp fails if atoms of the same type have different labels!
        # this only crops up as an error in the 4.0 versions
        # 'decimal_only' stops fractions that cannot be parsed when reading
        # updated positions
        if terse:
            # Don't output the bonds
            keywords = "opti noautobond decimal_only conj\n"
        else:
            keywords = "opti noautobond bond decimal_only conj\n"
        gin_file = [
            "# \n# Keywords:\n# \n",
            keywords,
            "# \n# Options:\n# \n",
            "# UFF optimisation by gulp\n"
            "name %s\n" % self.name,
            # using an rfo minimiser with a preconditioned hessian from the
            # BFGS minimiser seems to be the most efficient
            "switch rfo gnorm 0.3\n",
            "vectors\n"] + self.cell.to_vector_strings() + [
            " 1 1 1\n 1 1 1\n 1 1 1\n",  # constant pressure relaxation
            "cartesian\n"]

        all_ff_types = {}

        for at_idx, atom in enumerate(self.atoms):
            ff_type = atom.uff_type
            if not ff_type in all_ff_types:
                all_ff_types[ff_type] = atom.site
                # Sanity check in case types are not given in the input
                if ff_type not in UFF_FULL:
                    error("Atom %s has unknown type %s" % (atom, ff_type))
            if atom.is_fixed:
                fixed_flags = "0 0 0"
            else:
                fixed_flags = "1 1 1"
            gin_file.extend(["%-5s core " % all_ff_types[ff_type],
                             "%14.7f %14.7f %14.7f " % tuple(atom.pos),
                             "%f " % atom.charge,
                             "%f " % 1.0,  # occupancy
                             fixed_flags, "\n"])


        #identify all the individual uff species for the library
        gin_file.append("\nspecies\n")
        for ff_type, species in all_ff_types.items():
            gin_file.append("%-6s core %-6s\n" % (species, ff_type))

        gin_file.append("\n")
        for bond in sorted(self.bonds):
            bond_type = GULP_BOND_ORDERS[self.bonds[bond]]
            gin_file.append("connect %6i %6i %s\n" % (bond[0] + 1, bond[1] + 1, bond_type))

        gin_file.append("\nlibrary uff\n")

        gin_file.append("\nstepmx opt 0.05\n")

        # Restart file is for final structure
        gin_file.append("\ndump every %s.grs\n" % self.name)
        # optimization movie useful for debugging mostly
        if terse:
            # These tell gulp to be quiet, but we also stop the massive arc
            # file being generated
            gin_file.append("\nterse inout structure\n"
                            "terse inout potentials\n"
                            "terse inout derivatives\n"
                            "#output movie arc %s\n" % self.name)
        else:
            gin_file.append("\noutput movie arc %s\n" % self.name)

        return gin_file


    def to_egulp(self, typed_atoms=False):
        """Generate input files for Eugene's QEq code."""
        # bind cell locally for speed and convenience
        cell = self.cell.cell
        inv_cell = self.cell.inverse
        # format exactly as Eugene original script generates it
        geometry_file = ['%s\n' % self.name]
        geometry_file.extend(self.cell.to_vector_strings(fmt=' %15.12f'))
        geometry_file.append('%i\n' % self.natoms)

        atomic_numbers = self.atomic_numbers

        if typed_atoms:
            # Include custom typing, new types must be found by hand.
            # 800 N=O -- removed
            # 801 S=O
            # 802 S-O-H
            for atom_idx, atom in enumerate(self.atoms):
                if atom.uff_type == 'S_3+6':
                    for bond in self.bonds:
                        if atom_idx in bond:
                            other_idx = other_bond_index(bond, atom_idx)
                            if self.atoms[other_idx].uff_type == 'O_2':
                                atomic_numbers[other_idx] = 801
                            elif self.atoms[other_idx].uff_type == 'O_3':
                                atomic_numbers[other_idx] = 802
                                for another_bond in self.bonds:
                                    if other_idx in another_bond:
                                        another_idx = other_bond_index(another_bond, other_idx)
                                        if self.atoms[another_idx].uff_type == 'H_':
                                            atomic_numbers[another_idx] = 1001
# TODO(tdaff): delete code in future version; removed NO2 typing
#                elif atom.uff_type == 'N_R':
#                    this_bonds = []
#                    for bond in self.bonds:
#                        if atom_idx in bond:
#                            other_idx = other_bond_index(bond, atom_idx)
#                            if self.atoms[other_idx].uff_type == 'O_R':
#                                atomic_numbers[other_idx] = 800

        # atomic numbers should have been modified with exotic types by now
        geometry_file.extend([
            ('%6d ' % atomic_number) +
            ('%12.7f %12.7f  %12.7f' % tuple(atom.ipos(cell, inv_cell))) +
            ('%12.7f\n' % atom.charge)
            for atom, atomic_number in zip(self.atoms, atomic_numbers)
        ])

        return geometry_file


    def to_fastmc(self, options):
        """Return the FIELD and CONFIG needed for a fastmc run."""
        # CONFIG
        self.gen_supercell(options)
        supercell = self.gcmc_supercell
        levcfg = 0  # always
        imcon = self.cell.imcon
        natoms = len(self.atoms) * prod(supercell)
        config = ["%s\n" % self.name[:80],
                  "%10i%10i%10i\n" % (levcfg, imcon, natoms)]
        config.extend(self.cell.to_vector_strings(scale=supercell))
        for idx, atom in enumerate(self.supercell(supercell)):
            # idx+1 for 1 based indexes in CONFIG
            config.extend(["%-6s%10i\n" % (atom.type, idx + 1),
                           "%20.12f%20.12f%20.12f\n" % tuple(atom.pos)])

        # FIELD
        ntypes = len(self.guests) + 1
        field = ["%s\n" % self.name[:80],
                 "UNITS   kcal\n",
                 "molecular types %i\n" % ntypes]
        # Guests
        for guest in self.guests:
            field.extend(["&guest %s: %s\n" % (guest.name, guest.source),
                          "NUMMOLS %i\n" % 0,
                          "ATOMS %i\n" % len(guest.atoms)])
            for atom in guest.atoms:
                field.append(("%-6s %12.6f %12.6f" %
                              tuple([atom.type, atom.mass, atom.charge])) +
                             ("%12.6f %12.6f %12.6f\n" % tuple(atom.pos)))
            field.append("finish\n")
        # Framework
        field.extend(["Framework\n",
                      "NUMMOLS %i\n" % prod(supercell),
                      "ATOMS %i\n" % len(self.atoms)])
        if options.getbool('mc_zero_charges'):
            for atom in self.atoms:
                field.append("%-6s %12.6f %20.14f %6i %6i\n" %
                             (atom.type, atom.mass, 0.0, 1, 1))
        else:
            for atom in self.atoms:
                field.append("%-6s %12.6f %20.14f %6i %6i\n" %
                             (atom.type, atom.mass, atom.charge, 1, 1))
        field.append("finish\n")
        # VDW potentials
        atom_set = [atom.type for atom in self.atoms]
        for guest in self.guests:
            atom_set.extend(atom.type for atom in guest.atoms)
        atom_set = unique(atom_set)
        field.append("VDW %i\n" % ((len(atom_set) * (len(atom_set) + 1)) / 2))

        # modify local ff to deal with guests
        force_field = copy(UFF)
        for guest in self.guests:
            force_field.update(guest.potentials)

        for idxl in range(len(atom_set)):
            for idxr in range(idxl, len(atom_set)):
                left = atom_set[idxl]
                right = atom_set[idxr]
                try:
                    sigma, epsilon = lorentz_berthelot(force_field[left],
                                                       force_field[right])
                except KeyError:
                    # catch this if not in the UFF -> zero
                    warning("No potential defined for %s %s; defaulting to 0" %
                         (left, right))
                    sigma, epsilon = 0.0, 0.0
                field.append("%-6s %-6s lj %f %f\n" %
                             (left, right, epsilon, sigma))
        # EOF
        field.append("close\n")

        return config, field

    def to_gromacs(self):
        """Procedure:
        generate .gro atom positions
        generate .top lj topology
        generate .itp uff topology
        gererate .mdp parameters
        grompp -f mdp_file -c gro_file -p top_file -o output.tpr
        mdrun -s moffive -o traj-c moffive -nt 1
        """

        ##
        # .gro file
        ##
        # TITLE linw
        # number of atoms
        gro_file = ["%s\n" % self.name, " %i\n" % self.natoms]
        # TODO(tdaff): check residue in toplogy?
        residue = (1, self.name[:4].upper())
        for idx, atom in enumerate(self.atoms):
            # In the specification the atom position is given as
            # %8.3f in nm, but we probably would like more accuracy
            # but make sure there are 5 figures to left of decimal point
            pos_nm = tuple([x/10.0 for x in atom.pos])
            gro_file.append(("%5i%5s" % residue) +
                            ("%5s%5i" % (atom.uff_type, idx + 1)) +
                            ("%15.10f%15.10f%15.10f\n" % pos_nm))

        # cell is also in nm
        # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
        box = self.cell.cell/10.0
        # gromacs needs these condtiions. Might need to rotate cell sometimes?
        if not (box[0][1] == box[0][2] == box[1][2] == 0):
            error("Gromacs can't handle this cell orientation")
        gro_file.append("%f %f %f %f %f %f %f %f %f\n" % (
            box[0][0], box[1][1], box[2][2],
            box[0][1], box[0][2],
            box[1][0], box[1][2],
            box[2][0], box[2][1]))

        ##
        # .top file
        ##
        comb_rule = 1 # 1 = LORENTZ_BERTHELOT; 2 = GEOMETRIC
        top_file = [
            "; Topology file for %s\n" % self.name,
            "[ defaults ]\n",
            "; nbfunc      comb-rule      gen-pairs       fudgeLJ    fudgeQQ\n",
            "1             %i              yes             1.0        1.0\n\n" % (comb_rule)]

        # define the forcefield for UFF
        unique_types = set()
        top_file.append("[ atomtypes ]\n")
        top_file.append("; name1 name2   mass     charge  ptype   sigma   epsilon\n")
        for atom in self.atoms:
            if atom.uff_type not in unique_types:
                uff_type = atom.uff_type
                # sigma = x1 * 2^(-1/6)
                # CONVERT TO nm !!
                sigma = 0.1*UFF_FULL[uff_type][2]*(2**(-1.0/6.0))
                # epsilon = D1 in kcal
                epsilon = UFF_FULL[uff_type][3] * KCAL_TO_KJ
                top_file.append("%-6s   %-6s   %9.4f   %9.4f   A %12.7f %12.7f\n" % (uff_type, uff_type, atom.mass, 0.0, sigma, epsilon))
                unique_types.add(uff_type)

        top_file.extend(["\n#include <%s.itp>\n\n" % self.name,
                         "[ system ]\n",
                         "UFF optimisation of %s\n\n" % self.name,
                         "[ molecules ]\n",
                         "%s  1\n" % residue[1]])

        ##
        # .itp file
        ##
        # exclude 3 neighbours
        itp_file = ["[ moleculetype ]\n",
                    "; molname nrexcl\n"
                    "%s 3\n" % residue[1],
                    "[ atoms ]\n",
                    "; nr type  resnr    residue    atom     cgnr    charge       mass \n"]

        # atoms
        for idx, atom in enumerate(self.atoms):
            uff_type = atom.uff_type
            # charge group is different for each atom as gromacs has max of 32 in a group
            itp_file.append("%-6i  %-6s  %i  %-6s  %-6s   %d  %9.4f  %9.4f\n" % (idx+1, uff_type, residue[0], residue[1], uff_type, idx+1, atom.charge, atom.mass))

        # bonds
        itp_file.append(" [ bonds ]\n")
        itp_file.append("; ai aj funct b0 kb\n")

        unique_bonds = {}
        bonding_table = dict([(x, {}) for x in range(self.natoms)])
        #FIXME(tdaff) bonds will have length in v2
        for bond, bondorder in self.bonds.items():
            idx_a = bond[0]
            idx_b = bond[1]
            atom_a = self.atoms[idx_a]
            atom_b = self.atoms[idx_b]
            uff_a = atom_a.uff_type
            uff_b = atom_b.uff_type
            typed_bond = tuple(sorted([uff_a, uff_b]) + [bondorder])

            bonding_table[idx_a][idx_b] = typed_bond
            bonding_table[idx_b][idx_a] = typed_bond

            # have we already calculated the parameters
            if not typed_bond in unique_bonds:
                ri = UFF_FULL[uff_a][0]
                rj = UFF_FULL[uff_b][0]
                chiI = UFF_FULL[atom_a.uff_type][8]
                chiJ = UFF_FULL[atom_b.uff_type][8]
                rbo = -0.1332*(ri+rj)*log(bondorder)

                ren = ri*rj*(((sqrt(chiI) - sqrt(chiJ))**2)) / (chiI*ri + chiJ*rj)
                r0 = 0.1*(ri + rj + rbo - ren)

                # force constant
                # parameters Z1
                kb = (KCAL_TO_KJ * 664.12 * UFF_FULL[uff_a][5] * UFF_FULL[uff_b][5])/(r0**3)

                unique_bonds[typed_bond] = (r0, kb)

            params = unique_bonds[typed_bond]
            bond_func = 1  # gromacs harmonic
            # add 1 to bond as 1 indexed
            itp_file.append('%-5i %-5i %1i %11.4f %11.4f ; %-5s %-5s %.2f\n' % (bond[0]+1, bond[1]+1, bond_func, params[0], params[1], typed_bond[0], typed_bond[1], bondorder))

        # angles
        itp_file.append("\n [ angles ]\n")

        for idx_a in sorted(bonding_table):
            bonded_atoms = bonding_table[idx_a]
            for l_idx in bonded_atoms:
                for r_idx in bonded_atoms:
                    if l_idx == r_idx:
                        continue
                    # FIXME(tdaff) special cases for 5 or more coord
                    central_atom = self.atoms[idx_a]
                    l_atom = self.atoms[l_idx]
                    r_atom = self.atoms[r_idx]
                    coordination = len(bonded_atoms)
                    # FIXME(tdaff) over/under coordinated atoms?
                    theta0 = UFF_FULL[central_atom.uff_type][1]
                    cosT0 = cos(theta0*DEG2RAD)
                    sinT0 = sin(theta0*DEG2RAD)
                    c2 = 1.0 / (4.0 * sinT0 * sinT0)
                    c1 = -4.0 * c2 * cosT0
                    c0 = c2*(2.0*cosT0*cosT0 + 1.0)
                    zi = UFF_FULL[l_atom.uff_type][5]
                    zk = UFF_FULL[r_atom.uff_type][5]
                    bond_ab = tuple(sorted((l_idx, idx_a)))
                    bond_ab = tuple(sorted([l_atom.uff_type, central_atom.uff_type]) + [self.bonds[bond_ab]])
                    bond_bc = tuple(sorted((r_idx, idx_a)))
                    bond_bc = tuple(sorted([r_atom.uff_type, central_atom.uff_type]) + [self.bonds[bond_bc]])
                    rab = unique_bonds[bond_ab][0]
                    rbc = unique_bonds[bond_bc][0]
                    rac = sqrt(rab*rab + rbc*rbc - 2.0 * rab*rbc*cosT0)

                    ka = (644.12 * KCAL_TO_KJ) * (zi * zk / (rac**5.0))
                    ka *= (3.0*rab*rbc*(1.0 - cosT0*cosT0) - rac*rac*cosT0)

                    thetamin = (pi-arccos(c1/(4.0*c2))/DEG2RAD)
                    kappa = ka * (16.0*c2*c2 - c1*c1) / (4.0* c2)
                    itp_file.append("%-5i %-5i %-5i %3i %f %f ; %-5s %-5s %-5s" % (l_idx + 1, idx_a + 1, r_idx + 1, 1, thetamin, kappa, l_atom.uff_type, central_atom.uff_type, r_atom.uff_type))

        #with open("moffive.top", 'w') as tempfile:
        #    tempfile.writelines(top_file)
        #with open("moffive.gro", 'w') as tempfile:
        #    tempfile.writelines(gro_file)
        #with open("moffive.itp", 'w') as tempfile:
        #    tempfile.writelines(itp_file)


    def to_cssr(self, cartesian=False, no_atom_id=False):
        """
        Return a Cerius2 cssr file with coordinates and cell as a list of
        strings. Set no_atom_id to produce labels to work with Zeo++.

        """

        space_group = (1, "P 1")
        opt = 1
        cell = self.cell
        # spacers between all format strings ensures that there will always
        # be whitespace between components.
        cssr = [" "*38, "%8.3f %7.3f %7.3f\n" % cell.params[:3],
                " "*21, "%8.3f %7.3f %7.3f" % cell.params[3:], " "*4,
                "SPGR = %3i %-11s" % space_group, "OPT = %i\n" % opt,
                "%4i %3i %s\n" % (self.natoms, cartesian, self.name),
                "Structure file generated by faps\n"]

        for at_idx, atom in enumerate(self.atoms):
            # The default of TypeID does not work for Zeo++ as it
            # treats each label as a different type
            if no_atom_id:
                atom_name = "%s" % (atom.type)
            else:
                atom_name = "%s%i" % (atom.type, at_idx+1)
            # Zeo++ needs fractional coordinates
            if cartesian:
                string_pos = "%9.5f %9.5f %9.5f " % tuple(atom.pos)
            else:
                string_pos = "%9.5f %9.5f %9.5f " % tuple(atom.ifpos(cell.inverse))
            cssr.append("%4i %-4s  %s" % (at_idx+1, atom_name, string_pos) +
                        "   0"*8 + " %7.3f\n" % atom.charge)

        return cssr

    def to_zeoplusplus(self):
        """
        Return a tuple containing a cssr file, radii file and mass file.

        """

        radii = ["%-7s %-f\n" % (atom, UFF[atom][0]/2.0)
                 for atom in unique(self.types)]
        masses = ["%-7s %-f\n" % (atom, WEIGHT[atom])
                  for atom in unique(self.types)]

        return (self.to_cssr(no_atom_id=True), radii, masses)

    def fastmc_postproc(self, filepath, tp_point, options):
        """Update structure properties from gcmc OUTPUT."""
        startdir = os.getcwd()
        os.chdir(filepath)

        # Since Pete changed output strip indentation and blank lines
        filetemp = open('OUTPUT')
        output = strip_blanks(filetemp.readlines())
        filetemp.close()

        # Keep track of supercell so we can get unit cell values
        supercell_mult = prod(self.gcmc_supercell)
        # Still positional as we need multiple values simultaneously
        # and very old versions changed wording of heat of adsorption
        # and enthalpy of guest
        # TODO(tdaff, r2.0): deprecate reading older fastmc files
        # and put Cv in the guest definition
        for line in output[::-1]:
            if "+/-" in line:
                # This is version 1.3 of fastmc
                debug("NEW OUTPUT")
                line_offset = 5
                read_cv = True
                # In future this should be assumed to exist
                for guest in self.guests:
                    if not hasattr(guest, 'c_v'):
                        guest.c_v = {}
                break
        else:
            # +/- not found, assume old style output
            debug("OLD OUTPUT")
            line_offset = 0
            read_cv = False
        for idx, line in enumerate(output):
            # Assume that block will always start like this
            if 'final stats' in line:
                idx += line_offset
                guest_id = int(line.split()[4]) - 1
                self.guests[guest_id].uptake[tp_point] = (
                    float(output[idx + 1].split()[-1]),
                    float(output[idx + 2].split()[-1]),
                    supercell_mult)
                # This will sometimes be NaN
                self.guests[guest_id].hoa[tp_point] = (
                    float(output[idx + 3].split()[-1]),
                    float(output[idx + 4].split()[-1]))
                if read_cv:
                    self.guests[guest_id].c_v[tp_point] = (
                    float(output[idx + 5].split()[-1]),
                    float(output[idx + 6].split()[-1]))
            elif 'total accepted steps' in line:
                counted_steps = int(line.split()[-1])
                if counted_steps < 10000:
                    warning("Number of accepted GCMC steps is very low; "
                            "only %i counted!" % counted_steps)

        fold = options.getbool('fold')
        find_maxima = options.getbool('find_maxima')
        prob_plot = options.getbool('mc_probability_plot')
        folded = False
        if prob_plot and (fold or find_maxima):
            folded = self.fold_and_maxima(fold, find_maxima, tp_point)

        if folded and not options.getbool('fastmc_keep_unfolded_cubes'):
            debug("Removing unfolded cube files")
            cubes = glob.glob("prob_guest??_prob_??.cube")
            remove_files(cubes)

        #TODO(tdaff): absl

        unneeded_files = options.gettuple('fastmc_delete_files')
        remove_files(unneeded_files)
        keep_files = options.gettuple('fastmc_compress_files')
        compress_files(keep_files)

        os.chdir(startdir)


    def fold_and_maxima(self, fold=True, find_maxima=True, tp_point=None):
        """Determine the positions of maxima and produce an xyz xyz file."""
        from cube import Cube
        folded = False
        if fold:
            fold = self.gcmc_supercell
        else:
            fold = None
        for guest_idx, guest in enumerate(self.guests):
            guest_locations = {}
            for site_idx, sites in enumerate(guest.probability):
                guest_cube = Cube("prob_guest%02i_prob_%02i.cube" %
                                  (guest_idx+1, site_idx+1), fold=fold)
                if fold is not None:
                    debug("Folded cube file: %s" % guest_cube.folded_name)
                    guest_cube.write_cube()
                    folded = True
                if find_maxima:
                    guest_locations[sites] = guest_cube.maxima()
            if guest_locations:
                if tp_point:
                    # We can keep them for later too, must create dict
                    # since it might not exist in old calculations
                    if not hasattr(guest, 'guest_locations'):
                        guest.guest_locations = {tp_point: guest_locations}
                    else:
                        guest.guest_locations[tp_point] = guest_locations
                maxima = []
                for sites in guest_locations:
                    atom_name = name_from_types(sites, guest)
                    for atom, magnitude in guest_locations[sites]:
                        maxima.append((magnitude, "%-6s" % atom_name +
                                      "%10.6f %10.6f %10.6f " % tuple(atom) +
                                      "%10.6f\n" % magnitude))
                locations = open("%s-%s.xyz" % (self.name, guest.ident), 'w')
                locations.write(" %i\nBinding sites at %r\n" %
                                (len(maxima), tp_point))
                locations.writelines([x[1] for x in sorted(maxima, reverse=True)])
                locations.close()

        return folded


    def remove_duplicates(self, tolerance=0.02):
        """Find overlapping atoms and remove them."""
        uniq_atoms = []
        found_atoms = []
        for atom in self.atoms:
            for uniq_atom in uniq_atoms:
                if atom.type != uniq_atom.type:
                    continue
                elif min_distance(atom, uniq_atom) < tolerance:
                    break
            # else excutes when not found here
            else:
                uniq_atoms.append(atom)
        debug("Found %i unique atoms in %i" % (len(uniq_atoms), self.natoms))
        self.atoms = uniq_atoms

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

    def void_volume(self):
        """Estimate the void volume based on VdW radii."""
        initial_resolution = 0.2
        params = self.cell.params
        cell = self.cell.cell
        inv_cell = np.linalg.inv(cell.T)

        grid_size = [int(ceil(x/initial_resolution)) for x in params[:3]]
        print(grid_size)
        grid_resolution = [params[0]/grid_size[0],
                           params[1]/grid_size[1],
                           params[2]/grid_size[2]]
        print(grid_resolution)
        grid = np.zeros(grid_size, dtype=bool)
        atoms = [(atom.ipos(cell), inv_cell.tolist(),
                  atom.ifpos(inv_cell),
                  atom.vdw_radius) for atom in self.atoms]

        for x_idx in range(grid_size[0]):
            print(x_idx)
            for y_idx in range(grid_size[1]):
                print(y_idx)
                print(grid.sum())
                for z_idx in range(grid_size[2]):
                    grid_pos = [x_idx*grid_resolution[0],
                                y_idx*grid_resolution[1],
                                z_idx*grid_resolution[2]]
                    grid_fpos = [float(x_idx)/grid_size[0],
                                 float(y_idx)/grid_size[1],
                                 float(z_idx)/grid_size[2]]
                    for pos, fpos, radius in atoms:
                        dist = min_dist(pos, fpos, grid_pos, grid_fpos, cell)
                        if dist < radius:
                            grid[x_idx, y_idx, z_idx] = 1
                            break

        print(grid.sum())

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
    # TODO(tdaff): properties: density, surface area, dft_energy, absorbance

