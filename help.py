#!/usr/bin/env python

"""
PyTurds -- high throughput strucutre adsorption property analysis

As a script, run analysis on a structure. Provides classes and methods
for adapting the simulation or only doing select parts:

"""

import sys
import pickle
from optparse import OptionParser
from numpy import pi, cos, sin, sqrt, arccos
from numpy import array
from config import Options

DEG2RAD = pi/180

class Simulation(object):
    """
    A single property calculation for one structure.

    """
    # TODO(tdaff): automate the whole thing unless told otherwise
    def __init__(self, options):
        self.options = options
        self.structure = Structure("moffy")

    def run_vasp(self):
        pass


class Structure(object):
    """
    The current state of the structure; updates as the calculation proceeds.

    All simulation methods are defined for the structure as their input file
    formats etc.

    """
    # FIXME: no symmetry right now, eh?
    def __init__(self, name):
        self.name = name
        self.remark = False
        self.cell = Cell()
        self.atoms = []
        self.types = []
        self.esp = None
        # testing pdb reader for now
        self.from_pdb()

    def from_pdb(self, filename="MOF-5.pdb"):
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

    def _update_types(self):
        """Regenrate the list of atom types."""
        # FIXME: better ways of dealing with this; will come with mols/symmetry
        self.types = [atom.type for atom in self.atoms]

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


class Cell(object):
    """
    Crystollagraphic cell representations and interconversion methods.

    To ensure that cell and params are consistent, each format should use a
    setter method [ick] that calls the private interconversion methods.

    """
    def __init__(self):
        self.cell = array([[1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]])
        self.params = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)


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
        """[Super]cell vectors for dl_poly CONFIG."""
        return ["%20.12f%20.12f%20.12f\n" % tuple(scale*self.cell[0]),
                "%20.12f%20.12f%20.12f\n" % tuple(scale*self.cell[1]),
                "%20.12f%20.12f%20.12f\n" % tuple(scale*self.cell[2])]

    def to_vasp(self, scale=1):
        """[Super]cell vectors for VASP POSCAR."""
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



if __name__ == '__main__':
    global_options = Options()
    my_simulation = Simulation(global_options)
    #print(pickle.dumps(my_simulation))
