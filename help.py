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
from numpy import pi, cos, sin, sqrt, arccos
from numpy import array
from config import Options

DEG2RAD = pi/180
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
        for at_type in set(self.structure.types):
            potcar_src = os.path.join(self.options.potcar_dir, at_type, "POTCAR")
            shutil.copyfileobj(open(potcar_src), filetemp)
        filetemp.close()

        # TODO(tdaff): wooki specific at the moment
        vaspargs = ["vaspsubmit-beta", job_name, "%i" % nproc]
        submit = subprocess.Popen(vaspargs, stdout=subprocess.PIPE)
        for line in submit.stdout.readlines():
            if "wooki" in line:
                jobid = line.split(".")[0]
                print jobid
            else:
                print("Job failed?")



    def run_repeat(self):
        pass

    def run_cpmd(self):
        pass

    def run_fastmc(self):
        pass


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
    def __init__(self, name):
        """Just instance an empty structure initially"""
        self.name = name
        self.cell = Cell()
        self.atoms = []
        self.types = []
        self.esp = None
        # testing pdb reader for now
        self.from_pdb(self.name + '.pdb')

    def from_pdb(self, filename='MOF-5.pdb'):
        """Read an initial structure from a pdb file"""
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

        self._update_types()
        print self.cell.to_dl_poly(scale=2)
        print self.atoms
        sys.stdout.writelines(self.to_vasp())

    def from_cif(self, filename="structure.cif"):
        """Genereate structure from a .cif file"""
        raise NotImplementedError

    def charges_from_repeat(self, filename):
        charges = []
        filetemp = open(filename)
        for line in filetemp:
            if line.startswith(" Charge"):
                line = line.split()
                charges.append((int(line[1]), int(line[4]), float(line[6])))
            if "Error" in line:
                if float(line.split()[-1]) > 0.6:
                    print("Error in repeat charges is very high!")
        filetemp.close()
        # TODO: update structure

    def to_vasp(self, optim_h=True):
        """Return a vasp5 poscar as a list of lines."""
        poscar = ["%s\n" % self.name[:80],
                  " 1.0\n"]
        poscar.extend(self.cell.to_vasp())
        poscar.append("".join("%5s" % x for x in set(self.types)) + "\n")
        poscar.append("".join("%6i" % self.types.count(x) for x in set(self.types)) + "\n")
        if optim_h:
            poscar.extend(["Selective dynamics\n", "Cartesian\n"])
            for type in set(self.types):
                for atom in self.atoms:
                    if atom.type == type and type == "H":
                        poscar.append("%20.16f%20.16f%20.16f   T   T   T\n" % tuple(atom.pos))
                    elif  atom.type == type:
                        poscar.append("%20.16f%20.16f%20.16f   F   F   F\n" % tuple(atom.pos))
        else:
            poscar.append("Cartesian\n")
            for type in set(self.types):
                for atom in self.atoms:
                    if atom.type == type:
                        poscar.append("%20.16f%20.16f%20.16f\n" % tuple(atom.pos))
        return poscar

    def to_cpmd(self, optim_h=True):
        """Return a cpmd input file as a list of lines."""
        pass

    def to_fastmc(self):
        """Return the FIELD and CONFIG needed for a fastmc run"""
        pass

    def _update_types(self):
        """Regenrate the list of atom types."""
        # FIXME: better ways of dealing with this; will come with mols/symmetry
        self.types = [atom.type for atom in self.atoms]


class Cell(object):
    """
    Crystollagraphic cell representations and interconversion methods.

    To ensure that cell and params are consistent, each format should use a
    setter method [ick] that calls the private interconversion methods.

    """

    def __init__(self):
        """Default to a 1A box"""
        self.cell = array([[1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]])
        self.params = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

    def from_params(self, params):
        """Set the params and update the cell representation."""
        self.params = params
        self._mkcell()

    def from_pdb(self, line):
        """Extract cell from CRYST1 line in a pdb."""
        # TODO: space groups?
        self.params = tuple(float(x) for x in line.split()[1:7])
        self._mkcell()

    def _mkcell(self):
        """Update the cell representation to match the parameters."""
        a_mag, b_mag, c_mag = self.params[:3]
        alpha, beta, gamma = [x*DEG2RAD for x in self.params[3:]]
        a_vec = array([a_mag, 0.0, 0.0])
        b_vec = array([b_mag*cos(gamma), b_mag*sin(gamma), 0.0])
        c_x = c_mag*cos(beta)
        c_y = c_mag*(cos(alpha)-cos(gamma)*cos(beta))/sin(gamma)
        c_vec = array([c_x, c_y, (c_mag**2 - c_x**2 - c_y **2)**0.5])
        self.cell = array([a_vec, b_vec, c_vec])

    def _mkparam(self):
        """Update the parameters to match the cell."""
        cell_a = sqrt(sum(x**2 for x in self.cell[0:10:3])) # Sum of x
        cell_b = sqrt(sum(x**2 for x in self.cell[1:10:3]))
        cell_c = sqrt(sum(x**2 for x in self.cell[2:10:3]))
        alpha = arccos(sum(self.cell[1 + 3*j]*self.cell[2 + 3*j]
                           for j in range(3))/(cell_b*cell_c))*180/pi
        beta = arccos(sum(self.cell[0 + 3*j]*self.cell[2 + 3*j]
                          for j in range(3))/(cell_a*cell_c))*180/pi
        gamma = arccos(sum(self.cell[0 + 3*j]*self.cell[1 + 3*j]
                           for j in range(3))/(cell_a*cell_b))*180/pi
        self.params = (cell_a, cell_b, cell_c, alpha, beta, gamma)

    def to_dl_poly(self, scale=1):
        """[Super]cell vectors for DL_POLY CONFIG."""
        return ["%20.12f%20.12f%20.12f\n" % tuple(scale*self.cell[0]),
                "%20.12f%20.12f%20.12f\n" % tuple(scale*self.cell[1]),
                "%20.12f%20.12f%20.12f\n" % tuple(scale*self.cell[2])]

    def to_vasp(self, scale=1):
        """[Super]cell vectors for VASP POSCAR."""
        # Vasp usually has 16 dp but we get rounding errors eg cubic :(
        return ["%23.14f%22.14f%22.14f\n" % tuple(scale*self.cell[0]),
                "%23.14f%22.14f%22.14f\n" % tuple(scale*self.cell[1]),
                "%23.14f%22.14f%22.14f\n" % tuple(scale*self.cell[2])]

    def to_cpmd(self, scale=1):
        """[Super]cell vectors for CPMD input."""
        scale = scale/BOHR2ANG
        return ["%23.16f%22.16f%22.16f\n" % tuple(scale*self.cell[0]),
                "%23.16f%22.16f%22.16f\n" % tuple(scale*self.cell[1]),
                "%23.16f%22.16f%22.16f\n" % tuple(scale*self.cell[2])]


class Atom(object):
    """Base atom type, minimally requires type and position."""

    def __init__(self, at_type=False, pos=False):
        self.type = at_type
        self.pos = pos
        self.charge = 0
        self.idx = False
        self.site = None

    def __str__(self):
        return "%s %f %f %f" % tuple([self.type] + list(self.pos))

    def __repr__(self):
        return "Atom(%r,%r)" % (self.type, self.pos)

    def from_pdb(self, line):
        """Parse the ATOM line from a pdb file."""
        # Better to do this fixed width fields rather than splitting?
        self.idx = int(line[6:11])
        self.site = line[12:16].strip()
        self.molecule = int(line[22:26])
        at_pos = float(line[30:38]), float(line[38:46]), float(line[47:54])
        self.pos = at_pos
        self.type = line[76:78].strip()

    def translate(self, vec):
        """Move the atom by the given vector."""
        self.pos = [x + y for x, y in zip(self.pos, vec)]


def mk_repeat(cube_name='REPEAT_ESP.cube', symmetry=False):
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
    repeatsubmit()

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
