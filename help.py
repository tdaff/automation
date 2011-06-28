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
from copy import copy
from numpy import pi, cos, sin, sqrt, arccos, prod
from numpy import array, identity
from numpy import dot
from config import Options
from elements import WEIGHT, UFF

DEG2RAD = pi / 180.0
BOHR2ANG = 0.52917720859


class PyNiss(object):
    """
    PyNiss -- Negotiation of Intermediate System States

    A single property calculation for one structure.

    """
    # TODO(tdaff): automate the whole thing unless told otherwise
    def __init__(self, options):
        self.options = options
        self.structure = Structure(options.job_name)
        self.state = {'h_opt': 0,
                      'repeat': 0,
                      'gcmc': 0}

    def job_dispatcher(self):
        if self.options.interactive:
            import code
            console = code.InteractiveConsole(locals())
            console.interact()
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
        if self.options.h_optimize:
            print self.options.repeat_exe

    def run_vasp(self, nproc=16):
        """Make inputs and run vasp job."""
        job_name = self.options.job_name

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
            potcar_src = os.path.join(self.options.potcar_dir, at_type,
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
            else:
                print("Job failed?")

    def esp2cube(self):
        """Make the cube for repeat input."""
        job_name = self.options.job_name
        if self.options.esp_calc == 'vasp':
            os.chdir(job_name + ".restart_DIR")
            esp2cube_args = shlex.split(self.options.vasp2cube)
            submit = subprocess.Popen(esp2cube_args)
            # TODO(tdaff): maybe background this?
            submit.wait()
            # TODO(tdaff): leave the cube name as job-name..
            shutil.move(job_name + '.cube', self.options.cwd)
            os.chdir(self.options.cwd)

    def run_repeat(self):
        """Submit the repeat calc to the queue."""
        job_name = self.options.job_name
        mk_repeat(cube_name=job_name + '.cube')
        repeat_args = ['repeatsubmit', job_name + '.cube']
        submit = subprocess.Popen(repeat_args, stdout=subprocess.PIPE)
        for line in submit.stdout.readlines():
            if "wooki" in line:
                jobid = line.split(".")[0]
                print jobid
            else:
                print("Job failed?")

    def run_cpmd(self):
        pass

    def run_fastmc(self):
        confy, feeld = self.structure.to_fastmc(supercell=(2, 1, 1))


class Structure(object):
    """
    The current state of the structure; updates as the calculation proceeds.

    All simulation methods are defined for the structure as their input file
    formats etc.

    """
    # Methods are grouped:
    # * Structure parsers
    # * Output file parsers to update structure
    # * Input file generation
    # * internal manipulation methods
    # FIXME: no symmetry right now, eh?
    # TODO: dft energy?
    def __init__(self, name):
        """Just instance an empty structure initially."""
        self.name = name
        self.cell = Cell()
        self.atoms = []
        self.types = []
        self.esp = None
        self.dft_energy = 0.0
        self.guests = [Guest()]
        # testing pdb reader for now
        self.from_pdb(self.name + '.pdb')

    def from_pdb(self, filename='MOF-5.pdb'):
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

    def charges_from_repeat(self, filename):
        """Parse charges and update structure."""
        charges = []
        filetemp = open(filename)
        for line in filetemp:
            if line.startswith(" Charge"):
                line = line.split()
                charges.append((int(line[1]), int(line[4]), float(line[6])))
            if "Error" in line:
                if float(line.split()[-1]) > 0.6:
                    print("Error in repeat charges is very high - check cube!")
        filetemp.close()
        # TODO: update structure
        # TODO(tdaff): no symmetry here yet!
        for atom, charge in zip(self.atoms, charges):
            atom.charge = charge[2]

    def to_vasp(self, optim_h=True):
        """Return a vasp5 poscar as a list of lines."""
        ordered_types = unique(self.types)
        poscar = ["%s\n" % self.name[:80],
                  " 1.0\n"]
        poscar.extend(self.cell.to_vasp())
        poscar.append("".join("%5s" % x for x in ordered_types) + "\n")
        poscar.append("".join("%6i" % self.types.count(x)
                              for x in ordered_types) + "\n")
        if optim_h:
            poscar.extend(["Selective dynamics\n", "Cartesian\n"])
            for at_type in ordered_types:
                for atom in self.atoms:
                    if atom.type == at_type and at_type == "H":
                        poscar.append("%20.16f%20.16f%20.16f   T   T   T\n" %
                                      tuple(atom.pos))
                    elif  atom.type == at_type:
                        poscar.append("%20.16f%20.16f%20.16f   F   F   F\n" %
                                      tuple(atom.pos))
        else:
            poscar.append("Cartesian\n")
            for at_type in ordered_types:
                for atom in self.atoms:
                    if atom.type == at_type:
                        poscar.append("%20.16f%20.16f%20.16f\n" %
                                      tuple(atom.pos))
        return poscar

    def to_cpmd(self, optim_h=True):
        """Return a cpmd input file as a list of lines."""
        pass

    def to_fastmc(self, supercell=(1, 1, 1)):
        """Return the FIELD and CONFIG needed for a fastmc run"""
        # CONFIG
        # TODO(tdaff): guess imcon
        levcfg = 0
        imcon = 3
        natoms = len(self.atoms) * prod(supercell)
        atom_set = []
        config = [
            "%s\n" % self.name[:80],
            "%10i%10i%10i\n" % (levcfg, imcon, natoms)]
        config.extend(self.cell.to_vector_string(scale=supercell))
        for idx, atom in enumerate(self.supercell(supercell)):
            config.extend(["%-6s%10i\n" % (atom.type, idx),
                           "%20.12f%20.12f%20.12f\n" % tuple(atom.pos)])

        # FIELD
# TODO(tdaff): ntypes = nguests + nummols
        ntypes = 1 + len(self.guests)
        field = [
            "%s\n" % self.name[:80],
            "UNITS   kcal\n",
            "molecular types %i\n" % ntypes]
        # Guests
        for guest in self.guests:
            field.extend([
                "&guest %s\n" % guest.name,
                "NUMMOLS %i\n" % 0,
                "ATOMS %i\n" % len(guest.atoms)
            ])
            for atom in guest.atoms:
                field.append(("%-6s %12.6f %12.6f" %
                              tuple([atom.type, atom.mass, atom.charge])) +
                             ("%12.6f %12.6f %12.6f\n" % atom.pos))
            field.append("finish\n")
        # Framework
        field.extend([
            "Framework\n",
            "NUMMOLS %i\n" % prod(supercell),
            "ATOMS %i\n" % len(self.atoms)])
        for atom in self.atoms:
            field.append("%-6s %12.6f %20.14f %6i %6i\n" %
                         (atom.type, atom.mass, atom.charge, 1, 1))
        field.append("finish\n")
        # VDW potentials
        atom_set = self.types
        for guest in self.guests:
            atom_set.extend(guest.types)
        atom_set = unique(atom_set)  # TODO(tdaff) + guest types
        field.append("VDW %i\n" % ((len(atom_set) * (len(atom_set) + 1)) / 2))
        for idxl in range(len(atom_set)):
            for idxr in range(idxl, len(atom_set)):
                field.append(len_jones(atom_set[idxl], atom_set[idxr]))

        sys.stdout.writelines(field)

        return config, field

    def from_vasp(self, filename='CONTCAR', update=True):
        """Read a structure from a vasp [POS,CONT]CAR file."""
        #TODO(tdaff): difference between initial and update?
        filetemp = open(filename)
        contcar = filetemp.readlines()
        filetemp.close()
        atom_list = []
        scale = float(contcar[1])
        self.cell.from_vasp(contcar[2:4], scale)
        if contcar[5].split()[0].isalpha():
            # vasp 5 with atom names
#            self.types = []
            poscar_types = [x for x in contcar[5].split()]
#            for at_count, at_type in zip(atom_counts, poscar_types):
#                self.types.extend([at_type]*at_count)
            del contcar[5]
        else:
            #TODO atom ids when not in poscar?
            pass
        atom_counts = [int(x) for x in contcar[5].split()]
        if contcar[7].strip()[0].lower() in "s":
            # 's'elective dynamics line
            del contcar[7]
        # mcell converts frac -> cart and scales
        if contcar[7].strip()[0].lower() in "ck":
            mcell = identity(3) * scale
        else:
            mcell = self.cell
        if update:
            for atom in self.atoms:
                atom.from_vasp(at_line, cell=mcell)
        else:
            for at_type, at_line in zip(self.types, contcar[8:]):
                this_atom = Atom()
                this_atom.from_vasp(at_line, at_type, mcell)
                atom_list.append(this_atom)

        self.atoms = atom_list
        self._update_types()

    def supercell(self, dimensions):
        """Generate the atoms for a supercell"""
        if len(dimensions) == 1:
            dimensions = (dimensions, dimensions, dimensions)
        for x_super in range(dimensions[0]):
            for y_super in range(dimensions[1]):
                for z_super in range(dimensions[2]):
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
        self.cell = cell
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
        cell_a = sqrt(sum(x**2 for x in self.cell[0:10:3]))  # Sum of x
        cell_b = sqrt(sum(x**2 for x in self.cell[1:10:3]))
        cell_c = sqrt(sum(x**2 for x in self.cell[2:10:3]))
        alpha = arccos(sum(self.cell[1 + 3 * j] * self.cell[2 + 3 * j]
                           for j in range(3)) / (cell_b * cell_c)) * 180 / pi
        beta = arccos(sum(self.cell[0 + 3 * j] * self.cell[2 + 3 * j]
                          for j in range(3)) / (cell_a * cell_c)) * 180 / pi
        gamma = arccos(sum(self.cell[0 + 3 * j] * self.cell[1 + 3 * j]
                           for j in range(3)) / (cell_a * cell_b)) * 180 / pi
        self.params = (cell_a, cell_b, cell_c, alpha, beta, gamma)

    def to_vector_string(self, scale=1, bohr=False):
        """Generic [Super]cell vectors."""
        if bohr:
            scale = scale / BOHR2ANG
        return ["%20.12f%20.12f%20.12f\n" % tuple(scale * self.cell[0]),
                "%20.12f%20.12f%20.12f\n" % tuple(scale * self.cell[1]),
                "%20.12f%20.12f%20.12f\n" % tuple(scale * self.cell[2])]

    def to_dl_poly(self, scale=1):
        """[Super]cell vectors for DL_POLY CONFIG."""
        return ["%20.12f%20.12f%20.12f\n" % tuple(scale * self.cell[0]),
                "%20.12f%20.12f%20.12f\n" % tuple(scale * self.cell[1]),
                "%20.12f%20.12f%20.12f\n" % tuple(scale * self.cell[2])]

    def to_vasp(self, scale=1):
        """[Super]cell vectors for VASP POSCAR."""
        # Vasp usually has 16 dp but we get rounding errors eg cubic :(
        return ["%23.14f%22.14f%22.14f\n" % tuple(scale * self.cell[0]),
                "%23.14f%22.14f%22.14f\n" % tuple(scale * self.cell[1]),
                "%23.14f%22.14f%22.14f\n" % tuple(scale * self.cell[2])]

    def to_cpmd(self, scale=1):
        """[Super]cell vectors for CPMD input."""
        scale = scale / BOHR2ANG
        return ["%23.16f%22.16f%22.16f\n" % tuple(scale * self.cell[0]),
                "%23.16f%22.16f%22.16f\n" % tuple(scale * self.cell[1]),
                "%23.16f%22.16f%22.16f\n" % tuple(scale * self.cell[2])]


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
        """Set the atom data from vasp input"""
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


def mk_incar(job_name, h_opt=True):
    """Basic vasp INCAR; use defaults as much as possible."""
    incar = [
        "SYSTEM  = %s\n" % job_name,
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
    if h_opt:
        incar.extend([
            "IBRION  = 2\n",
            "NSW     = 300\n",
            "ISIF    = 2\n"])
    else:
        incar.extend([
            "IBRION  = 0\n",
            "NSW     = 0\n",
            "ISIF    = 0\n"])

    # TODO(tdaff)
        # "ENCUT = 520\n"
        # "#PREC    = high\n",
        # "ISPIN    = 2    ! Spin polarized plz\n",

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
        "steps                   40000000\n",
        "equilibration           4000000\n",
        "# jobcontrol\n",
        "cutoff          12.5 angstrom\n",
        "delr            1.0 angstrom\n",
        "ewald precision  1d-6\n",
        "finish\n"]
    return control


def len_jones(left, right):
    """Lorentz-Berthelot mixing rules for atom types"""
    sigma = (UFF[left][0] + UFF[right][0]) / 2.0
    epsilon = (UFF[left][1] * UFF[right][1])**0.5
    #TODO(tdaff): zero for zero?
    return "%-6s %-6s lj %f %f\n" % (left, right, epsilon, sigma)


if __name__ == '__main__':
    global_options = Options()
    # try to unpickle the job
    # fall back to starting a new simulation
    my_simulation = PyNiss(global_options)
    my_niss = open(global_options.job_name + ".niss", "wb")
    pickle.dump(my_simulation, my_niss)
    my_niss.close()

#    load_niss = open(global_options.job_name + ".niss")
#    my_simulation = pickle.load(load_niss)
#    load_niss.close()
#    print(my_simulation.structure.cell.to_dl_poly())
    my_simulation.job_dispatcher()
