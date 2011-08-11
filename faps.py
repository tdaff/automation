#!/usr/bin/env python

"""
PyFaps -- Fully Automated Pickles of Systems.

Strucutre adsorption property analysis for high throughput processing. When run
as a script, will automatically run complete analysis on a structure. Provides
classes and methods for adapting the simulation or only doing select parts.
Use the job.ini file to specify job-specific options; run each job in it's own
directory.

"""

__version__ = "$Revision$"

import code
import ConfigParser
import logging
import os
import pickle
import re
import shlex
import shutil
import subprocess
import sys
import textwrap
from copy import copy
from math import ceil

import numpy as np
from numpy import pi, cos, sin, sqrt, arccos, prod
from numpy import array, identity, dot, cross
from numpy.linalg import norm

from config import Options
from elements import WEIGHT, UFF, VASP_PSEUDO_PREF
from logo import LOGO

# Global constants
DEG2RAD = pi / 180.0
BOHR2ANG = 0.52917720859
EV2KCAL = 23.060542301389

# ID values for system state
NOT_RUN = 0
RUNNING = 1
FINISHED = 2
UPDATED = 3
SKIPPED = -1


class PyNiss(object):
    """
    PyNiss -- Negotiation of Intermediate System States

    A single property calculation for one structure. Instance with a set of
    options, then run the job_dispatcher() to begin the calculation. The
    calculation will pickle itself, or can be pickled at any time, by calling
    dump_state().

    """
    # TODO(tdaff): automate the whole thing unless told otherwise
    def __init__(self, options):
        """
        Instance an empty structure in the calculation; The dispatcher should
        be called to fill it up with data, as needed.

        """
        self.options = options
        self.structure = Structure(options.get('job_name'))
        self.state = {'init': (NOT_RUN, False),
                      'dft': (NOT_RUN, False),
                      'esp': (NOT_RUN, False),
                      'charges': (NOT_RUN, False),
                      'gcmc': (NOT_RUN, False)}

    def dump_state(self):
        """Write the .niss file holding the current system state."""
        job_name = self.options.get('job_name')
        info("Writing state file, %s.niss." % job_name)
        os.chdir(self.options.get('job_dir'))
        my_niss = open(job_name + ".niss", "wb")
        pickle.dump(self, my_niss)
        my_niss.close()

    def re_init(self, new_options):
        """Re initialize simulation (with updated options)."""
        if new_options.getbool('update_opts'):
            info("Using new options.")
            self.options = new_options
        else:
            # Just update command line stuff
            info("Using old options with new command line arguments.")
            self.options.args = new_options.args
            self.options.options = new_options.options
            self.options.cmdopts = new_options.cmdopts
        self.status(initial=True)

    def job_dispatcher(self):
        """
        Run parts explicity specified on the command line or do the next step
        in an automated run. Drop to interactive mode, if requested.

        """

        if 'status' in self.options.args:
            self.status(initial=True)

        if self.options.getbool('interactive'):
            console = code.InteractiveConsole(locals())
            console.interact(
                banner="""See manual for instructions on interactive use.""")

        if self.options.getbool('import'):
            info("Importing results from a previous simulation")
            self.import_old()
            self.dump_state()

        if self.state['init'][0] == NOT_RUN:
            info("Reading in structure")
            # No structure, should get one
            self.structure.from_file(
                self.options.get('job_name'),
                self.options.get('initial_structure_format'))
            self.state['init'] = (UPDATED, False)
            self.dump_state()

        # TODO(tdaff): does dft/optim always generate ESP?
        if self.state['dft'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_dft'):
                info("Skipping DFT step completely")
                warn("Job might fail if you need the ESP")
                self.state['dft'] = (SKIPPED, False)
            elif self.state['dft'][0] == RUNNING:
                new_state = jobcheck(self.state['dft'][1])
                if not new_state or new_state == 'C':
                    info("Queue reports DFT step has finished")
                    # Finished running update positions
                    self.structure.update_pos(self.options.get('dft_code'))
                    self.state['dft'] = (UPDATED, False)
                    self.dump_state()
                else:
                    # Still running
                    info("DFT still in progress")
                    terminate(0)

        if self.state['dft'][0] == NOT_RUN or 'dft' in self.options.args:
            self.run_dft()
            #info("Running optimizaton/dft step")
            self.dump_state()
            terminate(0)

        if self.state['charges'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_charges'):
                info("Skipping charge calculation")
                self.state['charges'] = (SKIPPED, False)
            elif self.state['charges'][0] == RUNNING:
                new_state = jobcheck(self.state['charges'][1])
                if not new_state or new_state == 'C':
                    info("Queue reports charge calculation has finished")
                    self.structure.update_charges(
                        self.options.get('charge_method'))
                    self.state['charges'] = (UPDATED, False)
                    self.dump_state()
                else:
                    info("Charge calculation still running")
                    terminate(0)

        if self.state['charges'][0] == NOT_RUN or 'charges' in self.options.args:
            self.run_charges()
            info("Running charge calculation")
            self.dump_state()
            terminate(0)

        if self.state['gcmc'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_gcmc'):
                info("Skipping GCMC simulation")
                self.state['gcmc'] = (UPDATED, False)
            elif self.state['gcmc'][0] == RUNNING:
                new_state = jobcheck(self.state['gcmc'][1])
                if not new_state or new_state == 'C':
                    info("Queue reports GCMC simulation has finished")
                    # Read in GCMC data
                    self.state['gcmc'] = (UPDATED, False)
                    self.dump_state()
                else:
                    info("GCMC still running")
                    terminate(0)

        if self.state['gcmc'][0] == NOT_RUN or 'gcmc' in self.options.args:
            self.run_fastmc()
            info("Running gcmc step")
            self.dump_state()
            terminate(0)
        else:
            # Everything finished
            info("GCMC run has finished")

    def status(self, initial=False):
        """Print the current status to the terminal."""
        valid_states = {NOT_RUN: 'Not run',
                        RUNNING: 'Running',
                        FINISHED: 'Finished',
                        UPDATED: 'Processed',
                        SKIPPED: 'Skipped'}

        if initial:
            info("Previous system state reported from .niss file "
                 "(running jobs may have already finished):")
        else:
            info("Current system status:")
        for step, state in self.state.iteritems():
            if state[0] is RUNNING:
                info(" * State of %s: Running, jobid: %s" % (step, state[1]))
            else:
                info(" * State of %s: %s" % (step, valid_states[state[0]]))

    def import_old(self):
        """Try and import any data from previous stopped simulation."""
        job_name = self.options.get('job_name')
        job_dir = self.options.get('job_dir')
        try:
            self.structure.from_file(
                job_name,
                self.options.get('initial_structure_format'))
            self.state['init'] = (UPDATED, False)
        except IOError:
            info("No initial structure found to import")
        try:
            self.structure.update_pos(self.options.get('dft_code'))
            self.state['dft'] = (UPDATED, False)
        except IOError:
            info("No optimized structure found to import")
        try:
            self.structure.update_charges(self.options.get('charge_method'))
            self.state['charges'] = (UPDATED, False)
        except IOError:
            info("No charges found to import")
        # Reset directory at end
        os.chdir(job_dir)

    def run_dft(self):
        """Select correct method for running the dft/optim."""
        dft_code = self.options.get('dft_code')
        info("Running a %s calculation" % dft_code)
        if dft_code == 'vasp':
            self.run_vasp()
        elif dft_code == 'cpmd':
            err("CPMD calculation not yet implemented")
        else:
            err("Unknown dft method")

    def run_charges(self):
        """Select correct charge processing methods."""
        chg_method = self.options.get('charge_method')
        info("Calculating charges with %s" % chg_method)
        if chg_method == 'repeat':
            # Get ESP
            self.esp2cube()
            self.run_repeat()
        else:
            err("Unknown charge calculation method: %s" % chg_method)

    def run_vasp(self):
        """Make inputs and run vasp job."""
        job_name = self.options.get('job_name')
        nproc = self.options.getint('vasp_ncpu')
        # Keep things tidy in a subdirectory
        dft_code = self.options.get('dft_code')
        vasp_dir = os.path.join(self.options.get('job_dir'),
                                'faps_%s_%s' % (job_name, dft_code))
        mkdirs(vasp_dir)
        os.chdir(vasp_dir)
        debug("Running in %s" % vasp_dir)
        info("Running on %i nodes" % nproc)

        filetemp = open("POSCAR", "wb")
        filetemp.writelines(self.structure.to_vasp(self.options))
        filetemp.close()

        filetemp = open("INCAR", "wb")
        filetemp.writelines(mk_incar(self.options))
        filetemp.close()

        filetemp = open("KPOINTS", "wb")
        filetemp.writelines(mk_kpoints(self.options.gettuple('kpoints')))
        filetemp.close()

        potcar_types = unique([atom.type for atom in self.structure.atoms])
        filetemp = open("POTCAR", "wb")
        for at_type in potcar_types:
            # Try and get the preferred POTCARS
            # TODO(tdaff): update these with custom pseudos
            debug("Using %s pseudopotential for %s" %
                 (VASP_PSEUDO_PREF.get(at_type, at_type), at_type))
            potcar_src = os.path.join(self.options.get('potcar_dir'),
                                      VASP_PSEUDO_PREF.get(at_type, at_type),
                                      "POTCAR")
            shutil.copyfileobj(open(potcar_src), filetemp)
        filetemp.close()

        # different names for input files on wooki
        for vasp_file in ['POSCAR', 'INCAR', 'KPOINTS', 'POTCAR']:
            shutil.copy(vasp_file, job_name + '.' + vasp_file.lower())

        if self.options.getbool('no_submit'):
            info("Vasp input files generated; skipping job submission")
            self.state['dft'] = (SKIPPED, False)
        else:
            # FIXME(tdaff): wooki specific at the moment
            vasp_args = ["vaspsubmit-beta", job_name, "%i" % nproc]
            submit = subprocess.Popen(vasp_args, stdout=subprocess.PIPE)
            for line in submit.stdout.readlines():
                if "wooki" in line:
                    jobid = line.split(".")[0]
                    info("Running VASP job in queue. Jobid: %s" % jobid)
                    self.state['dft'] = (RUNNING, jobid)
                    break
            else:
                warn("Job failed?")
        # Tidy up at the end
        os.chdir(self.options.get('job_dir'))

    def esp2cube(self):
        """Make the cube for repeat input."""
        job_name = self.options.get('job_name')
        esp_src = self.options.get('esp_src')
        src_dir = os.path.join(self.options.get('job_dir'),
                               'faps_%s_%s' % (job_name, esp_src))
        os.chdir(src_dir)
        if esp_src == 'vasp':
            os.chdir(job_name + ".restart_DIR")
            esp2cube_args = shlex.split(self.options.get('vasp2cube'))
            info("Converting esp to cube, this might take a minute...")
            submit = subprocess.Popen(esp2cube_args)
            submit.wait()
            # TODO(tdaff): leave the cube name as job-name..
            # Move it to the repeat directory
            repeat_dir = os.path.join(self.options.get('job_dir'),
                                      'faps_%s_repeat' % job_name)
            mkdirs(repeat_dir)
            move_and_overwrite(job_name + '.cube', repeat_dir)
            os.chdir(self.options.get('job_dir'))

    def run_repeat(self):
        """Submit the repeat calc to the queue."""
        job_name = self.options.get('job_name')
        charge_code = self.options.get('charge_method')
        repeat_dir = os.path.join(self.options.get('job_dir'),
                                  'faps_%s_%s' % (job_name, charge_code))
        mkdirs(repeat_dir)
        os.chdir(repeat_dir)

        mk_repeat(cube_name=job_name + '.cube')
        if self.options.getbool('no_submit'):
            info("REPEAT input files generated; skipping job submission")
            self.state['charges'] = (SKIPPED, False)
        else:
            repeat_args = ['repeatsubmit', job_name + '.cube']
            submit = subprocess.Popen(repeat_args, stdout=subprocess.PIPE)
            jobid = None
            for line in submit.stdout.readlines():
                if "wooki" in line:
                    jobid = line.split(".")[0]
                    info("Running REPEAT calculation in queue: Jobid %s"
                         % jobid)
                    self.state['charges'] = (RUNNING, jobid)
                    break
            else:
                warn("Job failed?")
        os.chdir(self.options.get('job_dir'))

    def run_cpmd(self):
        pass

    def run_fastmc(self):
        """Submit a fastmc job to the queue."""
        job_name = self.options.get('job_name')
        mc_code = self.options.get('mc_code')

        gcmc_dir = os.path.join(self.options.get('job_dir'),
                                'faps_%s_%s' % (job_name, mc_code))
        mkdirs(gcmc_dir)
        os.chdir(gcmc_dir)

        # Set the guests before generating the files
        # Load here as options may change in each run
        guests = self.options.get('guests')
        # Try and deal with any list of guests
        self.structure.guests = [Guest(x) for x in
                                 re.split('[\s,\(\)\[\]]*', guests) if x]
        config, field = self.structure.to_fastmc(self.options)

        filetemp = open("CONFIG", "wb")
        filetemp.writelines(config)
        filetemp.close()

        filetemp = open("FIELD", "wb")
        filetemp.writelines(field)
        filetemp.close()

        filetemp = open("CONTROL", "wb")
        filetemp.writelines(mk_gcmc_control(self.options,
                                            self.structure.guests))
        filetemp.close()

        if self.options.getbool('no_submit'):
            info("FastMC input files generated; skipping job submission")
            self.state['gcmc'] = (SKIPPED, False)
        else:
            fastmc_args = ['fastmcsubmit', job_name]
            submit = subprocess.Popen(fastmc_args, stdout=subprocess.PIPE)
            jobid = None
            for line in submit.stdout.readlines():
                if "wooki" in line:
                    jobid = line.split(".")[0]
                    info("Running FastMC in queue: Jobid %s" % jobid)
                    self.state['gcmc'] = (RUNNING, jobid)
                    break
            else:
                warn("Job submission failed?")
        os.chdir(self.options.get('job_dir'))


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
    # FIXME: no symmetry right now, eh?
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

    def from_file(self, basename, filetype):
        """Select the correct file parser."""
        filetype = filetype.lstrip('.')
        if filetype.lower() in ['pdb']:
            self.from_pdb(basename + '.' + filetype)
        elif filetype.lower() in ['vasp', 'poscar', 'contcar']:
            self.from_vasp()
        elif filetype.lower() in ['cif']:
            self.from_cif(basename + '.' + filetype)
        else:
            err("Unknown filetype %s" % filetype)

    def update_pos(self, opt_code):
        """Select the method for updating atomic positions."""
        opt_path = os.path.join('faps_%s_%s' % (self.name, opt_code))
        info("Updating positions from %s" % opt_code)
        if opt_code == 'vasp':
            self.from_vasp(os.path.join(opt_path, self.name + '.contcar'),
                           update=True)
        elif opt_code == 'cpmd':
            self.from_cpmd(update=True)

    def update_charges(self, charge_method):
        """Select the method for updating charges."""
        charge_path = os.path.join('faps_%s_%s' % (self.name, charge_method))
        if charge_method == 'repeat':
            info("Updating charges from repeat")
            self.charges_from_repeat(
                os.path.join(charge_path, self.name + '.esp_fit.out'))

    def from_pdb(self, filename):
        """Read an initial structure from a pdb file."""
        info("Reading positions from pdb file: %s" % filename)
        filetemp = open(filename)
        pdb_file = filetemp.readlines()
        filetemp.close()
        # Build a local list before setting attribute
        newatoms = []
        for line in pdb_file:
            if line.lower().startswith('cryst1'):
                self.cell.from_pdb(line)
            elif line.lower().startswith('atom'):
                newatom = Atom()
                newatom.from_pdb(line)
                newatoms.append(newatom)
            elif line.lower().startswith('hetatm'):
                newatom = Atom()
                newatom.from_pdb(line)
                newatoms.append(newatom)

        self.atoms = newatoms
        self.order_by_types()

    def from_cif(self, filename):
        """Genereate structure from a .cif file."""
        info("Reading positions from cif file: %s" % filename)
        filetemp = open(filename)
        cif_file = filetemp.readlines()
        filetemp.close()
        cif_file = without_blanks(cif_file)
        params = [None]*6
        atoms = []
        symmetry = []
        loops = []
        idx = 0
        while idx < len(cif_file):
            line = cif_file[idx].lower()
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
            elif 'loop_' in line:
                # loops for _atom_site and _symmetry
                heads = []
                body = []
                while '_' in line:
                    heads.extend(line.split())
                    idx += 1
                    line = cif_file[idx]
                while idx < len(cif_file) and '_' not in line:
                    # shlex keeps 'quoted items' as one
                    body.extend(shlex.split(line))
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
            err("No cell or incomplete cell found in cif file")

        # TODO: symmetry
        # parse loop contents
        for heads, body in loops:
            if '_atom_site_fract_x' in heads:
                while body:
                    atoms.append(dict(zip(heads, body)))
                    body = body[len(heads):]
            if '_symmetry_equiv_pos_as_xyz' in heads:
                while body:
                    sym_dict = dict(zip(heads, body))
                    symmetry.append(Symmetry(sym_dict['_symmetry_equiv_pos_as_xyz']))
                    body = body[len(heads):]

        newatoms = []
        for atom in atoms:
            for sym_op in symmetry:
                newatom = Atom()
                newatom.from_cif(atom, self.cell.cell, sym_op)
                newatoms.append(newatom)

        self.atoms = newatoms
        self.order_by_types()

    def from_vasp(self, filename='CONTCAR', update=True):
        """Read a structure from a vasp [POS,CONT]CAR file."""
        #TODO(tdaff): difference between initial and update?
        info("Reading positions from vasp file: %s" % filename)
        filetemp = open(filename)
        contcar = filetemp.readlines()
        filetemp.close()
        atom_list = []
        scale = float(contcar[1])
        self.cell.from_vasp(contcar[2:5], scale)
        if contcar[5].split()[0].isalpha():
            # vasp 5 with atom names
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
        else:
            # FIXME(tdaff): not working
            for at_type, at_line in zip(self.types, contcar[7:7+natoms]):
                this_atom = Atom()
                this_atom.from_vasp(at_line, at_type, mcell)
                atom_list.append(this_atom)
            self.atoms = atom_list

    def charges_from_repeat(self, filename):
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
                    warn("Error in repeat charges is very high - check cube!")
        filetemp.close()
        # TODO(tdaff): no symmetry here yet!
        for atom, charge in zip(self.atoms, charges):
            atom.charge = charge[2]

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
        for atom in self.atoms:
            if atom.type == "H":
                poscar.append("%20.16f%20.16f%20.16f" % tuple(atom.pos) +
                              "%4s%4s%4s\n" % (fix_h, fix_h, fix_h))
            else:
                poscar.append("%20.16f%20.16f%20.16f" % tuple(atom.pos) +
                              "%4s%4s%4s\n" % (fix_all, fix_all, fix_all))
        return poscar

    def to_cpmd(self, options):
        """Return a cpmd input file as a list of lines."""
        raise NotImplementedError

    def to_fastmc(self, options):
        """Return the FIELD and CONFIG needed for a fastmc run"""
        # TODO(tdaff): initialize guests
        # CONFIG
        config_supercell = options.gettuple('mc_supercell')
        config_cutoff = options.getfloat('mc_cutoff')
        if config_cutoff < 12:
            warn("Simulation is using a very small cutoff! I hope you "
                 "know, what you are doing!")
        minimum_supercell = self.cell.minimum_supercell(config_cutoff)
        self.gcmc_supercell = minimum_supercell
        supercell = tuple(max(i, j)
                          for i, j in zip(config_supercell, minimum_supercell))
        info("%s supercell requested in config" % str(config_supercell))
        info("%s minimum supercell for a %.1f cutoff" %
             (str(minimum_supercell), config_cutoff))
        info("Constructing %s supercell for gcmc." % str(supercell))
        levcfg = 0  # always
        imcon = self.cell.imcon
        natoms = len(self.atoms) * prod(supercell)
        config = ["%s\n" % self.name[:80],
                  "%10i%10i%10i\n" % (levcfg, imcon, natoms)]
        config.extend(self.cell.to_vector_strings(scale=supercell))
        for idx, atom in enumerate(self.supercell(supercell)):
            # idx+1 for 1 based indexes in CONFIG
            config.extend(["%-6s%10i\n" % (atom.type, idx+1),
                           "%20.12f%20.12f%20.12f\n" % tuple(atom.pos)])

        # FIELD
        # TODO(tdaff): ntypes = nguests + nummols
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
                             ("%12.6f %12.6f %12.6f\n" % atom.pos))
            field.append("finish\n")
        # Framework
        field.extend(["Framework\n",
                      "NUMMOLS %i\n" % prod(supercell),
                      "ATOMS %i\n" % len(self.atoms)])
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
        # TODO(tdaff): extra ones from options?
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
                    warn("No potential defined for %s %s; defaulting to 0" %
                         (left, right))
                    sigma, epsilon = 0.0, 0.0
                field.append("%-6s %-6s lj %f %f\n" %
                             (left, right, epsilon, sigma))
        # EOF
        field.append("close\n")

        return config, field

    def supercell(self, scale):
        """Iterate over all the atoms of supercell."""
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
        # TODO(tdaff): different for molecules

    @property
    def types(self):
        """Ordered list of atom types."""
        return [atom.type for atom in self.atoms]

    @property
    def weight(self):
        """Unit cell weight."""
        return sum([atom.mass for atom in self.atoms])

    @property
    def volume(self):
        """Unite cell volume."""
        return self.cell.volume

    def get_supercell(self):
        """Supercell used for gcmc."""
        if 'supercell' in self.properties:
            return self.properties['supercell']
        else:
            return None

    def set_supercell(self, value):
        """Set the supercell property for the structure."""
        self.properties['supercell'] = value

    gcmc_supercell = property(get_supercell, set_supercell)
    # TODO(tdaff): properties: volume, supercell, density, surface area
    # dft_energy, absorbance

class Cell(object):
    """
    Crystollagraphic cell representations and interconversion methods.

    Setter methods can be defined for different file types, however
    .cell and .params will be self-consistent if set directly.

    """

    def __init__(self):
        """Default to a 1A cubic box."""
        self._cell = array([[1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0]])
        self._params = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

    def from_pdb(self, line):
        """Extract cell from CRYST1 line in a pdb."""
        # TODO: space groups?
        # Must use fixed widths as -ve numbers do not leave gaps to .split()
        self.params = (float(line[6:15]),
                       float(line[15:24]),
                       float(line[24:33]),
                       float(line[33:40]),
                       float(line[40:47]),
                       float(line[47:54]))

    def from_vasp(self, lines, scale=1.0):
        """Extract cell from a POSCAR cell representation."""
        self.cell = array([[float(x) * scale for x in lines[0].split()],
                           [float(x) * scale for x in lines[1].split()],
                           [float(x) * scale for x in lines[2].split()]])

    def to_vector_strings(self, scale=1, bohr=False, fmt="%20.12f"):
        """Generic [Super]cell vectors in Angstrom as a list of strings."""
        out_format = 3 * fmt + "\n"
        if isinstance(scale, int):
            scale = [scale, scale, scale]
            # else assume an iterable
        if bohr:
            scale = [x / BOHR2ANG for x in scale]
        return [out_format % tuple(scale[0] * self.cell[0]),
                out_format % tuple(scale[1] * self.cell[1]),
                out_format % tuple(scale[2] * self.cell[2])]

    def minimum_supercell(self, cutoff):
        """Calculate the smallest supercell with a half-cell width cutoff."""

        a_cross_b = cross(self.cell[0], self.cell[1])
        b_cross_c = cross(self.cell[1], self.cell[2])
        c_cross_a = cross(self.cell[2], self.cell[0])

        volume = dot(self.cell[0], b_cross_c)

        widths = [volume/norm(b_cross_c),
                  volume/norm(c_cross_a),
                  volume/norm(a_cross_b)]

        return tuple(int(ceil(2*cutoff/x)) for x in widths)

    @property
    def imcon(self):
        """Guess cell shape and return DL_POLY imcon key."""
        keys = {'none': 0,
                'cubic': 1,
                'orthorhombic': 2,
                'parallelepiped': 3}
        if np.all(self.cell == 0):
            return keys['none']
        elif np.all(self.params[3:] == 90):
            if self.params[0] == self.params[1] == self.params[2]:
                return keys['cubic']
            else:
                return keys['orthorhombic']
        else:
            return keys['parallelepiped']

    @property
    def volume(self):
        """Calculate cell volume a.bxc."""
        b_cross_c = cross(self.cell[1], self.cell[2])
        return dot(self.cell[0], b_cross_c)

    def get_cell(self):
        """Get the three vector cell representation."""
        return self._cell
    def set_cell(self, value):
        """Set cell and params from the cell representation."""
        self._cell = value
        self.__mkparam()
    # Property so that params are updated simultaneously
    cell = property(get_cell, set_cell)

    def get_params(self):
        """Get the six parameter cell representation."""
        return self._params
    def set_params(self, value):
        """Set cell and params from the cell parameters."""
        self._params = value
        self.__mkcell()
    # Property so that cell is updated simultaneously
    params = property(get_params, set_params)

    # Implementation details -- directly access the private _{cell|param}
    # attributes; please don't break.
    def __mkcell(self):
        """Update the cell representation to match the parameters."""
        a_mag, b_mag, c_mag = self.params[:3]
        alpha, beta, gamma = [x * DEG2RAD for x in self.params[3:]]
        a_vec = array([a_mag, 0.0, 0.0])
        b_vec = array([b_mag * cos(gamma), b_mag * sin(gamma), 0.0])
        c_x = c_mag * cos(beta)
        c_y = c_mag * (cos(alpha) - cos(gamma) * cos(beta)) / sin(gamma)
        c_vec = array([c_x, c_y, (c_mag**2 - c_x**2 - c_y**2)**0.5])
        self._cell = array([a_vec, b_vec, c_vec])

    def __mkparam(self):
        """Update the parameters to match the cell."""
        cell_a = sqrt(sum(x**2 for x in self.cell[0]))
        cell_b = sqrt(sum(x**2 for x in self.cell[1]))
        cell_c = sqrt(sum(x**2 for x in self.cell[2]))
        alpha = arccos(sum(self.cell[1, :] * self.cell[2, :]) /
                       (cell_b * cell_c)) * 180 / pi
        beta = arccos(sum(self.cell[0, :] * self.cell[2, :]) /
                      (cell_a * cell_c)) * 180 / pi
        gamma = arccos(sum(self.cell[0, :] * self.cell[1, :]) /
                       (cell_a * cell_b)) * 180 / pi
        self._params = (cell_a, cell_b, cell_c, alpha, beta, gamma)


class Atom(object):
    """Base atom type."""

    def __init__(self, at_type=False, pos=False, **kwargs):
        """Accept arbritary kwargs as attributes."""
        self.type = at_type
        self.pos = pos
        self.charge = 0.0
        self.idx = False
        self.site = None
        self.mass = 0.0
        self.molecule = None
        # Sets anything else specified as an attribute
        for key, val in kwargs.iteritems():
            setattr(self, key, val)

    def __str__(self):
        return "%s %f %f %f" % tuple([self.type] + list(self.pos))

    def __repr__(self):
        return "Atom(%r,%r)" % (self.type, self.pos)

    def from_cif(self, at_dict, cell, symmetry=None):
        """Extract an atom description from dictionary of cif items."""
        self.type = at_dict['_atom_site_type_symbol']
        self.site = at_dict['_atom_site_label']
        self.mass = WEIGHT[self.type]
        frac_pos = [ufloat(at_dict['_atom_site_fract_x']),
                    ufloat(at_dict['_atom_site_fract_y']),
                    ufloat(at_dict['_atom_site_fract_z'])]
        if symmetry is not None:
            frac_pos = symmetry.trans_frac(frac_pos)
        self.pos = dot(frac_pos, cell)

    def from_pdb(self, line):
        """Parse the ATOM line from a pdb file."""
        # pdb is defined with fixed width fields rather than splitting
        self.idx = int(line[6:11])
        self.site = line[12:16].strip()
        self.molecule = int(line[22:26])
        at_pos = float(line[30:38]), float(line[38:46]), float(line[47:54])
        self.pos = at_pos
        self.type = line[76:78].strip()
        self.mass = WEIGHT[self.type]

    def from_vasp(self, line, at_type=None, cell=identity(3)):
        """Set the atom data from vasp input. Only pass cell if fractional."""
        self.pos = dot([float(x) for x in line.split()[:3]], cell)
        if at_type is not None:
            self.type = at_type
            self.mass = WEIGHT[at_type]

    def translate(self, vec):
        """Move the atom by the given vector."""
        self.pos = [x + y for x, y in zip(self.pos, vec)]


class Guest(object):
    """Guest molecule and properties."""
    def __init__(self, ident):
        self.ident = ident
        self.name = "Unknown guest"
        self.potentials = {}
        self.probability = []
        self.raw_guest = {}
        self.atoms = []
        self.source = "Unknown source"
        self._find_guest()
        self._parse_guest()

    def _find_guest(self):
        """Look in guests.lib in submit directory and default."""
        # Need the different directories
        job_dir = os.getcwd()
        if __name__ != '__main__':
            script_dir = os.path.dirname(__file__)
        else:
            script_dir = os.path.abspath(sys.path[0])
        job_guests = ConfigParser.SafeConfigParser()
        lib_guests = ConfigParser.SafeConfigParser()
        # Try and find guest in guests.lib
        job_guests.read(os.path.join(job_dir, 'guests.lib'))
        lib_guests.read(os.path.join(script_dir, 'guests.lib'))
        if job_guests.has_section(self.ident):
            debug("%s found in job dir" % self.ident)
            self.raw_guest = job_guests.items(self.ident)
        elif lib_guests.has_section(self.ident):
            debug("%s found in library" % self.ident)
            self.raw_guest = lib_guests.items(self.ident)
        else:
            err("Guest not found: %s" % self.ident)

    def _parse_guest(self):
        """Set attributes according to the raw input."""
        for key, val in self.raw_guest:
            if key == 'atoms':
                # Only use non blank lines
                atoms = [x.strip() for x in val.splitlines() if x.strip()]
                for atom in atoms:
                    atom = atom.split()
                    self.atoms.append(Atom(
                        type=atom[0],
                        mass=float(atom[1]),
                        charge=float(atom[2]),
                        pos=tuple(float(x) for x in atom[3:6])))
            elif key == 'potentials':
                potens = [x.strip() for x in val.splitlines() if x.strip()]
                for poten in potens:
                    poten = poten.split()
                    self.potentials[poten[0]] = tuple(float(x)
                                                      for x in poten[1:])
            elif key == 'probability':
                prob = [x.strip() for x in val.splitlines() if x.strip()]
                self.probability = [tuple(int(y) for y in x.split())
                                    for x in prob]
            else:
                setattr(self, key, val)


class Symmetry(object):
    """Apply symmetry operations to atomic coordinates."""
    def __init__(self, blob):
        """Read the operation from the argument."""
        self.sym_ops = []
        self.blob = blob
        self.parse_blob()
        self.cell = None

    def parse_blob(self):
        """Interpret a symmetry line from a cif."""
        # convert integers to floats to avoid integer division
        self.sym_ops = [re.sub(r'([\d]+)', r'\1.0', x.strip())
                        for x in re.split(',', self.blob) if x.strip()]

    def trans_frac(self, pos):
        """Apply symmetry operation to the supplied position."""
        # TODO(tdaff): check for overlapping atoms?
        new_pos = [eval(sym_op.replace('x', str(pos[0]))
                        .replace('y', str(pos[1]))
                        .replace('z', str(pos[2]))) for sym_op in self.sym_ops]
        return new_pos



def mk_repeat(cube_name='REPEAT_ESP.cube', symmetry=False):
    """Standard REPEAT input file."""
    # TODO(tdaff): charged systems?
    if symmetry:
        symmetry_flag = 1
        # TODO(tdaff): connectivity.ff
    else:
        symmetry_flag = 0
    repeat_input = [
        "Input ESP file name in cube format\n",
        "%s\n" % cube_name,
        "Fit molecular(0) or periodic(1:default) system?\n",
        "1\n",
        "van der Waals scaling factor (default = 1.0)\n",
        "1.00000\n",
        "Apply RESP penalties?, no(0:default), yes(1)\n",
        "0\n",
        "Read cutoff radius? no(0), yes(1:default)\n",
        "1\n",
        "If flag above=1 provide R_cutoff next (in Bohrs)\n",
        "20.00000\n",
        "Apply symmetry restrain? no(0:default), yes(1)\n",
        "%i\n" % symmetry_flag,
        "Use Goddard-type restrains? no(0:default), yes(1)\n",
        "0\n",
        "If flag above=1 then provide weight next\n",
        "0.00000\n",
        "Enter total charge of the system\n",
        "0.00000\n"]

    filetemp = open('REPEAT_param.inp', 'w')
    filetemp.writelines(repeat_input)
    filetemp.close()


def mk_incar(options):
    """Basic vasp INCAR; use defaults as much as possible."""
    # We need these options
    job_name = options.get('job_name')
    spin = options.getbool('spin')
    optim_h = options.getbool('optim_h')
    optim_all = options.getbool('optim_all')
    optim_cell = options.getbool('optim_cell')
    dispersion = options.getbool('dispersion')

    incar = ["SYSTEM  = %s\n" % job_name,
             "ALGO    = Fast\n",
             "EDIFF   = 1E-5\n",
             "EDIFFG  = -0.02\n",
             "POTIM   = 0.4\n",
             "LREAL   = Auto\n",
             "LVTOT   = .TRUE.\n",
             "LVHAR   = .TRUE.\n",
             "ISMEAR  = 0\n",
             "SIGMA   = 0.05\n"]
    if optim_cell:
        # Positions will be fixed by selective dynamics
        info("Cell vectors will be optimized")
        incar.extend(["ENCUT = 520\n",
                      "IBRION  = 2\n",
                      "NSW     = 300\n",
                      "ISIF    = 3\n"])
    elif optim_all or optim_h:
        # Just move positions
        incar.extend(["IBRION  = 2\n",
                      "NSW     = 300\n",
                      "ISIF    = 2\n"])
    else:
        # Single point energy
        info("Single point calculation")
        incar.extend(["IBRION  = 0\n",
                      "NSW     = 0\n",
                      "ISIF    = 0\n"])
    if spin:
        info("Spin polarised calculation")
        incar.append("ISPIN   = 2\n")

    if dispersion:
        info("Dispersion correction will be used")
        incar.append("LVDW    = .TRUE.\n")

    return incar


def mk_kpoints(kpoints):
    """Defaults to gamma point only, or specified number."""
    if len(kpoints) != 3:
        err("kpoints specified incorectly; should be (i, i, i)")
    kpoints = [
        "Auto\n",
        "0\n",
        "Gamma\n",
        "%i %i %i\n" % tuple(kpoints),
        "0 0 0\n"]
    return kpoints


def mk_gcmc_control(options, guests):
    """Standard GCMC CONTROL file."""
    control = [
        "GCMC Run\n"
        "temperature  %f\n" % options.getfloat('mc_temperature'),
        "steps    %i\n" % options.getint('mc_prod_steps'),
        "equilibration    %i\n" % options.getint('mc_eq_steps'),
        "cutoff          %f angstrom\n" % options.getfloat('mc_cutoff'),
        "delr            1.0 angstrom\n",
        "ewald precision  1d-6\n",
        "numguests %i\n" % options.getint('mc_numguests_freq'),
        "history %i\n" % options.getint('mc_history_freq')]
    if options.getbool('mc_jobcontrol'):
        control.append("jobcontrol\n")
    else:
        control.append("# jobcontrol\n")
    # Guest stuff
    if options.getbool('mc_probability_plot'):
        # We just comment out the probabilty plot lines if unneeded
        probp = ""
    else:
        probp = "# "
    guest_count = 0
    for guest in guests:
        guest_count += 1
        control.append("&guest %i\n" % guest_count)
        control.append("  pressure (bar) %f\n" %
                       options.getfloat('mc_pressure'))
        control.append("%s  probability %i\n" % (probp, len(guest.probability)))
        for prob in guest.probability:
            control.append("%s  %i  " % (probp, len(prob)) +
                           "  ".join(["%i" % x for x in prob]) + "\n")
        control.append("&end\n")

    control.append("finish\n")
    return control


def unique(in_list):
    """Unique values in list ordered by first occurance"""
    uniq = []
    for item in in_list:
        if item not in uniq:
            uniq.append(item)
    return uniq


def lorentz_berthelot(left, right):
    """Lorentz-Berthelot mixing rules for (sigma, epsilon) tuples."""
    sigma = (left[0] + right[0]) / 2.0
    epsilon = (left[1] * right[1])**0.5
    if epsilon == 0:
        sigma = 0
    return sigma, epsilon


def jobcheck(jobid):
    """Get job status."""
    jobid = "%s" % jobid
    qstat = subprocess.Popen(['qstat', '%s' % jobid], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
    for line in qstat.stdout.readlines():
        debug(line)
        if "Unknown Job Id" in line:
            return False
        elif line.startswith(jobid):
            status = line[68:69]
            return status
    else:
        print("Failed to get job information.")


def debug(msg):
    """Send DEBUGging to the logging handlers."""
    msg = textwrap.wrap(msg)
    for line in msg:
        logging.debug(line)


def info(msg):
    """Send INFO to the logging handlers."""
    msg = textwrap.wrap(msg)
    for line in msg:
        logging.info(line)


def warn(msg):
    """Send WARNings to the logging handlers."""
    msg = textwrap.wrap(msg)
    for line in msg:
        logging.warning(line)


def err(msg):
    """Send ERRORs to the logging handlers."""
    msg = textwrap.wrap(msg)
    for line in msg:
        logging.error(line)
    # TODO(tdaff): should we quit here?


def mkdirs(directory):
    """Create a directory if it does not exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)


def terminate(exit_code=0):
    """Exit and announce if faps is terminating normally (default)."""
    if exit_code == 0:
        info("Faps terminated normally")
        raise SystemExit
    else:
        warn("Abnormal termination of faps; exit code %i" % exit_code)
        raise SystemExit(exit_code)


def move_and_overwrite(src, dest):
    """Move src to dest and overwrite if it is an existing file."""
    if os.path.exists(dest):
        if os.path.isdir(dest):
            dest_full = os.path.join(dest, os.path.basename(src))
            if os.path.exists(dest_full):
                if os.path.isfile(dest_full):
                    os.remove(dest_full)
                    shutil.move(src, dest)
                else:
                    raise OSError("Directory %s already exists" % dest_full)
            else:
                shutil.move(src, dest)
        elif os.path.isfile(dest):
            os.remove(dest)
            shutil.move(src, dest)
        else:
            raise OSError("%s is not a folder or file" % dest)
    else:
        shutil.move(src, dest)

def ufloat(text):
    """Convert string to float, ignoring the uncertainty part."""
    return float(re.sub('\(.*\)', '', text))


def without_blanks(lines):
    """Strip lines and remove blank lines."""
    return [line.strip() for line in lines if line.strip() != '']


def welcome():
    """Print any important messages."""
    print("FAPS version 0.0r%s" % __version__.strip('$Revision: '))
    print(LOGO)
    print("\nFaps is under heavy development and will break without notice.")
    print("\nThis version breaks backwards and forwards compatibility!")
    print("EXISTING JOBS WILL NOT WORK [automatically], sorry :(")
    print("\n * Jobs now run in subdirectories;")
    print("   You need to manually move some files to --import existing jobs!")
    print("\n * Internal structure of the .niss file has changed,")
    print("   this should be deleted before '--import'ing.")
    print("\n")


def main():
    """Do a standalone calculation when run as a script."""
    welcome()
    main_options = Options()
    info("Starting FAPS version 0.0r%s" % __version__.strip('$Revision: '))
    # try to unpickle the job or
    # fall back to starting a new simulation
    niss_name = main_options.get('job_name') + ".niss"
    if os.path.exists(niss_name):
        info("Existing simulation found: %s; loading..." % niss_name)
        load_niss = open(niss_name)
        my_simulation = pickle.load(load_niss)
        load_niss.close()
        my_simulation.re_init(main_options)
    else:
        info("Starting a new simulation...")
        my_simulation = PyNiss(main_options)

    # run requested jobs
    my_simulation.job_dispatcher()
    my_simulation.dump_state()
    terminate(0)


if __name__ == '__main__':
    main()
