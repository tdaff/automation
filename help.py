#!/usr/bin/env python

"""
PyFaps -- Fully Automated Pickles of Systems.

High throughput strucutre adsorption property analysis. When run as a script,
will automatically run complete analysis on a structure. Provides classes and
methods for adapting the simulation or only doing select parts:

RUN -- run or continue until job is finished
STEP -- run only the next step and stop (specify multiple times for more steps)
STATUS -- print some information about the state of the system

valid steps:
XXX optimise and ESP
XXX new postions
XXX repeat
XXX new charges
XXX gcmc
XXX postproc

"""

import sys
import pickle
import shutil
import os
import subprocess
import shlex
import time
import numpy as np
from copy import copy
from numpy import pi, cos, sin, sqrt, arccos, prod
from numpy import array, identity
from numpy import dot
from config import Options
from elements import WEIGHT, UFF

DEG2RAD = pi / 180.0
BOHR2ANG = 0.52917720859
EV2KCAL = 23.060542301389

NOT_RUN = 0
RUNNING = 1
FINISHED = 2


class PyNiss(object):
    """
    PyNiss -- Negotiation of Intermediate System States

    A single property calculation for one structure.

    """
    # TODO(tdaff): automate the whole thing unless told otherwise
    def __init__(self, options):
        """An empty structure; The dispatcher will fill it up with data."""
        self.options = options
        self.structure = Structure(options.get('job_name'))
        self.state = {'init': (NOT_RUN, False),
                      'dft': (NOT_RUN, False),
                      'repeat': (NOT_RUN, False),
                      'gcmc': (NOT_RUN, False)}

    def job_dispatcher(self):
        """
        Drop to interactive mode if requested. Run parts explicity specified
        or the next step in automated run.

        """

        if self.options.get('interactive'):
            import code
            console = code.InteractiveConsole(locals())
            console.interact()
        if self.state['init'][0] == NOT_RUN:
            # No structure, should get one
            self.structure.from_file(options.get('job_name'),
                                     options.get('initial_structure_type'))
        if "opt" in self.options.args:
            self.run_vasp()
        if "esp" in self.options.args:
            self.esp2cube()
        if "repeat" in self.options.args:
            self.run_repeat()
        if "gcmc" in self.options.args:
            self.run_fastmc()
        if self.state['gcmc'] == 2:
            # Everything finished
            print("GCMC run has finished")
        if self.options.get('h_optimize'):
            print self.options.get('repeat_exe')

    def status(self):
        """Print the current status to the terminal."""
        print(self.status)

    def run_vasp(self, nproc=16):
        """Make inputs and run vasp job."""
        job_name = self.options.get('job_name')

        filetemp = open(job_name + ".poscar", "wb")
        filetemp.writelines(self.structure.to_vasp())
        filetemp.close()

        filetemp = open(job_name + ".incar", "wb")
        filetemp.writelines(mk_incar(job_name))
        filetemp.close()

        filetemp = open(job_name + ".kpoints", "wb")
        filetemp.writelines(mk_kpoints())
        filetemp.close()

        filetemp = open(job_name + ".potcar", "wb")
        for at_type in unique(self.structure.types):
            potcar_src = os.path.join(self.options.get('potcar_dir'), at_type,
                                      "POTCAR")
            shutil.copyfileobj(open(potcar_src), filetemp)
        filetemp.close()

        # TODO(tdaff): wooki specific at the moment
        vasp_args = ["vaspsubmit-beta", job_name, "%i" % nproc]
        submit = subprocess.Popen(vasp_args, stdout=subprocess.PIPE)
        for line in submit.stdout.readlines():
            if "wooki" in line:
                jobid = line.split(".")[0]
                print jobid
                self.state['dft'] = (RUNNING, jobid)
                break
        else:
            print("Job failed?")

    def esp2cube(self):
        """Make the cube for repeat input."""
        job_name = self.options.get('job_name')
        if self.options.get('esp_calc') == 'vasp':
            os.chdir(job_name + ".restart_DIR")
            esp2cube_args = shlex.split(self.options.get('vasp2cube'))
            submit = subprocess.Popen(esp2cube_args)
            # TODO(tdaff): maybe background this?
            submit.wait()
            # TODO(tdaff): leave the cube name as job-name..
            shutil.move(job_name + '.cube', self.options.get('cwd'))
            os.chdir(self.options.get('cwd'))

    def run_repeat(self):
        """Submit the repeat calc to the queue."""
        job_name = self.options.get('job_name')
        mk_repeat(cube_name=job_name + '.cube')
        repeat_args = ['repeatsubmit', job_name + '.cube']
        submit = subprocess.Popen(repeat_args, stdout=subprocess.PIPE)
        jobid = None
        for line in submit.stdout.readlines():
            if "wooki" in line:
                jobid = line.split(".")[0]
                print jobid
#        if jobid:
                self.state['repeat'] = (RUNNING, jobid)
                break
        else:
            print("Job failed?")

    def run_cpmd(self):
        pass

    def run_fastmc(self):
        """Submit a fastmc job to the queue."""
        confy, feeld = self.structure.to_fastmc(supercell=(1, 2, 1))


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
        self.types = []  # TODO(tdaff): remove and just use sorted self.atoms?
        self.esp = None
        self.dft_energy = 0.0
        self.guests = [Guest()]
        # FIXME(tdaff): just testing pdb reader for now
        #self.from_pdb(self.name + '.pdb')
        #self.from_vasp(self.name + '.contcar')
#        self.charges_from_repeat(self.name + '.esp_fit.out')

    def from_file(self, basename, filetype):
        """Select the correct file parser."""
        filetype = filetype.lstrip('.')
        if filetype.lower() in ['pdb']:
            self.from_pdb(basename + '.' + filetype)
        elif filetype.lower() in ['vasp', 'poscar', 'contcar']:
            self.from_vasp()


    def from_pdb(self, filename):
        """Read an initial structure from a pdb file."""
        filetemp = open(filename)
        pdbfile = filetemp.readlines()
        filetemp.close()
        for line in pdbfile:
            if 'cryst1' in line.lower():
                # the cell the input
                self.cell.from_pdb(line)
            if 'atom' in line.lower():
                newatom = Atom()
                newatom.from_pdb(line)
                self.atoms.append(newatom)

        self.oreder_by_types()
        self._update_types()

    def from_cif(self, filename="structure.cif"):
        """Genereate structure from a .cif file."""
        raise NotImplementedError

    def from_vasp(self, filename='CONTCAR', update=True):
        """Read a structure from a vasp [POS,CONT]CAR file."""
        #TODO(tdaff): difference between initial and update?
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
            for at_type, at_line in zip(self.types, contcar[7:7+natoms]):
                this_atom = Atom()
                this_atom.from_vasp(at_line, at_type, mcell)
                atom_list.append(this_atom)
            self.atoms = atom_list
        self._update_types()

    def charges_from_repeat(self, filename):
        """Parse charges and update structure."""
        charges = []
        filetemp = open(filename)
        for line in filetemp:
            if line.startswith(" Charge"):
                line = line.split()
                # index, type, charge
                charges.append((int(line[1]), int(line[4]), float(line[6])))
            if "Error" in line:
                if float(line.split()[-1]) > 0.6:
                    print("Error in repeat charges is very high - check cube!")
        filetemp.close()
        # TODO(tdaff): no symmetry here yet!
        for atom, charge in zip(self.atoms, charges):
            atom.charge = charge[2]

    def to_vasp(self, optim_h=True, optim_all=False):
        """Return a vasp5 poscar as a list of lines."""
        types = [atom.type for atom in self.atoms]
        ordered_types = unique(types)
        poscar = ["%s\n" % self.name[:80],
                  " 1.0\n"]
        # Vasp does 16 dp but we get rounding errors eg cubic -> 14
        poscar.extend(self.cell.to_vector_strings(fmt="%23.14f"))
        poscar.append("".join("%5s" % x for x in ordered_types) + "\n")
        poscar.append("".join("%6i" % types.count(x)
                              for x in ordered_types) + "\n")
        if optim_all:
            poscar.extend(["Selective dynamics\n", "Cartesian\n"])
            fix_h = 'T'
            fix_all = 'T'
        elif optim_h:
            poscar.extend(["Selective dynamics\n", "Cartesian\n"])
            fix_h = 'T'
            fix_all = 'F'
        else:
            poscar.append("Cartesian\n")
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

    def to_cpmd(self, optim_h=True):
        """Return a cpmd input file as a list of lines."""
        raise NotImplementedError

    def to_fastmc(self, supercell=(1, 1, 1)):
        """Return the FIELD and CONFIG needed for a fastmc run"""
        # CONFIG
        levcfg = 0  # always
        imcon = self.cell.imcon()
        natoms = len(self.atoms) * prod(supercell)
        config = ["%s\n" % self.name[:80],
                  "%10i%10i%10i\n" % (levcfg, imcon, natoms)]
        config.extend(self.cell.to_vector_strings(scale=supercell))
        for idx, atom in enumerate(self.supercell(supercell)):
            config.extend(["%-6s%10i\n" % (atom.type, idx),
                           "%20.12f%20.12f%20.12f\n" % tuple(atom.pos)])

        sys.stdout.writelines(config)

        # FIELD
        # TODO(tdaff): ntypes = nguests + nummols
        ntypes = 1 + len(self.guests)
        field = ["%s\n" % self.name[:80],
                 "UNITS   kcal\n",
                 "molecular types %i\n" % ntypes]
        # Guests
        for guest in self.guests:
            field.extend(["&guest %s\n" % guest.name,
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
        for idxl in range(len(atom_set)):
            for idxr in range(idxl, len(atom_set)):
                field.append(len_jones(atom_set[idxl], atom_set[idxr]))

        sys.stdout.writelines(field)

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

    def _update_types(self):
        """Regenrate the list of atom types."""
        # FIXME: better ways of dealing with this; will come with mols/symmetry
        self.types = [atom.type for atom in self.atoms]

    def oreder_by_types(self):
        """Sort the atoms alphabetically and group them."""
        self.atoms.sort(key=lambda x: (x.type, x.site))
        # TODO(tdaff): different for molecules


class Cell(object):
    """
    Crystollagraphic cell representations and interconversion methods.

    To ensure that cell and params are consistent, each format should use a
    setter method [ick] that calls the private interconversion methods.

    """

    def __init__(self):
        """Default to a 1A cubic box."""
        self.cell = array([[1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]])
        self.params = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

    def from_params(self, params):
        """Set the params and update the cell representation."""
        self.params = params
        self._mkcell()

    def from_cell(self, cell):
        """Set the params and update the cell representation."""
        self.cell = array(cell)
        self._mkparam()

    def from_pdb(self, line):
        """Extract cell from CRYST1 line in a pdb."""
        # TODO: space groups?
        self.params = tuple(float(x) for x in line.split()[1:7])
        self._mkcell()

    def from_vasp(self, lines, scale=1.0):
        """Extract cell from a POSCAR cell representation."""
        self.cell = array([[float(x) * scale for x in lines[0].split()],
                           [float(x) * scale for x in lines[1].split()],
                           [float(x) * scale for x in lines[2].split()]])
        self._mkparam()

    def _mkcell(self):
        """Update the cell representation to match the parameters."""
        a_mag, b_mag, c_mag = self.params[:3]
        alpha, beta, gamma = [x * DEG2RAD for x in self.params[3:]]
        a_vec = array([a_mag, 0.0, 0.0])
        b_vec = array([b_mag * cos(gamma), b_mag * sin(gamma), 0.0])
        c_x = c_mag * cos(beta)
        c_y = c_mag * (cos(alpha) - cos(gamma) * cos(beta)) / sin(gamma)
        c_vec = array([c_x, c_y, (c_mag**2 - c_x**2 - c_y**2)**0.5])
        self.cell = array([a_vec, b_vec, c_vec])

    def _mkparam(self):
        """Update the parameters to match the cell."""
        cell_a = sqrt(sum(x**2 for x in self.cell[0]))
        cell_b = sqrt(sum(x**2 for x in self.cell[1]))
        cell_c = sqrt(sum(x**2 for x in self.cell[2]))
        alpha = arccos(sum(self.cell[:, 1] * self.cell[:, 2]) /
                       (cell_b * cell_c)) * 180 / pi
        beta = arccos(sum(self.cell[:, 0] * self.cell[:, 2]) /
                      (cell_a * cell_c)) * 180 / pi
        gamma = arccos(sum(self.cell[:, 0] * self.cell[:, 1]) /
                       (cell_a * cell_b)) * 180 / pi
        self.params = (cell_a, cell_b, cell_c, alpha, beta, gamma)

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

    def imcon(self):
        """Guess cell shape and return DL_POLY imcon key."""
        if np.all(self.cell == 0):
            # no PBC
            return 0
        elif np.all(self.params[3:] == 90):
            if self.params[0] == self.params[1] == self.params[2]:
                # cubic
                return 1
            else:
                # orthorhombic
                return 2
        else:
            # parallelepiped
            return 3


class Atom(object):
    """Base atom type."""

    def __init__(self, at_type=False, pos=False, charge=0.0):
        self.type = at_type
        self.pos = pos
        self.charge = charge
        self.idx = False
        self.site = None
        self.mass = 0.0
        self.molecule = None

    def __str__(self):
        return "%s %f %f %f" % tuple([self.type] + list(self.pos))

    def __repr__(self):
        return "Atom(%r,%r)" % (self.type, self.pos)

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
        """Set the atom data from vasp input. Only pass cell for fractionals."""
        self.pos = dot([float(x) for x in line.split()[:3]], cell)
        if at_type is not None:
            self.type = at_type
            self.mass = WEIGHT[at_type]

    def translate(self, vec):
        """Move the atom by the given vector."""
        self.pos = [x + y for x, y in zip(self.pos, vec)]


class Guest(object):
    """Guest molecule and properties."""
    def __init__(self):
        self.name = "Carbon Dioxide"
        self.atoms = [
            Atom("Cx", pos=(0.0000000, 0.000000, 0.000000), charge=0.65120),
            Atom("Ox", pos=(1.1605000, 0.000000, 0.000000), charge=-0.32560),
            Atom("Ox", pos=(-1.1605000, 0.000000, 0.000000), charge=-0.32560)]
        self.types = [atom.type for atom in self.atoms]


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


def mk_incar(job_name, pos_opt=True, full_opt=False, spin=False):
    """Basic vasp INCAR; use defaults as much as possible."""
    incar = ["SYSTEM  = %s\n" % job_name,
             "ALGO    = Fast\n",
             "EDIFF   = 1E-5\n",
             "EDIFFG  = -0.02\n",
             "POTIM   = 0.4\n",
             "LVDW    = .TRUE.\n",
             "NWRITE  = 0\n",
             "LREAL   = Auto\n",
             "LVTOT   = .TRUE.\n",
             "LVHAR   = .TRUE.\n",
             "ISMEAR  = 0\n",
             "SIGMA   = 0.05\n"]
    if full_opt:
        incar.extend(["ENCUT = 520\n",
                      "IBRION  = 2\n",
                      "NSW     = 300\n",
                      "ISIF    = 3\n"])
    elif pos_opt:
        incar.extend(["IBRION  = 2\n",
                      "NSW     = 300\n",
                      "ISIF    = 2\n"])
    else:
        incar.extend(["IBRION  = 0\n",
                      "NSW     = 0\n",
                      "ISIF    = 0\n"])
    if spin:
        incar.append("ISPIN   = 2\n")

    return incar


def mk_kpoints(num_kpt=1):
    """Defaults to gamma point only, or specified number."""
    kpoints = [
        "Auto\n",
        "0\n",
        "Gamma\n",
        "%i %i %i\n" % (num_kpt, num_kpt, num_kpt),
        "0 0 0\n"]
    return kpoints


def unique(in_list):
    """Set of unique values in list ordered by first occurance"""
    uniq = []
    for item in in_list:
        if item not in uniq:
            uniq.append(item)
    return uniq


def mk_gcmc_control(num_guests, pressure):
    """Standard GCMC CONTROL file."""
    control = [
        "temperature               273\n",
        "&guest 1\n",
        "  pressure  (bar)  1.0\n",
        "&end\n",
        "steps                   1000000\n",
        "equilibration           100000\n",
        "# jobcontrol\n",
        "cutoff          12.5 angstrom\n",
        "delr            1.0 angstrom\n",
        "ewald precision  1d-6\n",
        "finish\n"]
    return control


def len_jones(left, right):
    """Lorentz-Berthelot mixing rules for atom types"""
    # UFF is found in elements.py
    # TODO(tdaff): better way of defaults/custom
    sigma = (UFF[left][0] + UFF[right][0]) / 2.0
    epsilon = (UFF[left][1] * UFF[right][1])**0.5
    if epsilon == 0:
        sigma = 0
    return "%-6s %-6s lj %f %f\n" % (left, right, epsilon, sigma)


if __name__ == '__main__':
    global_options = Options()
    # try to unpickle the job or
    # fall back to starting a new simulation
    if os.path.exists(global_options.get('job_name') + ".niss"):
        print("Existing simulation found; loading...")
        load_niss = open(global_options.get('job_name') + ".niss")
        my_simulation = pickle.load(load_niss)
        load_niss.close()
        my_simulation.options = global_options
    else:
        print("Starting a new simulation...")
        my_simulation = PyNiss(global_options)

    # run requested jobs
    my_simulation.job_dispatcher()

    # dump the final system state
    my_niss = open(global_options.get('job_name') + ".niss", "wb")
    pickle.dump(my_simulation, my_niss)
    my_niss.close()
