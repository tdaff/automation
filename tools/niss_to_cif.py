#!/usr/bin/env python

"""
niss_to_cif.py

Extract the structure from running or finished faps simulataions. Will process
multiple jobs when given together on the commandline. This script will read in
the state of the simulation from the niss file and create a cif with the latest
optimised structure including the most recent charges.

Usage: niss_to_cif.py JOB_NAME1 [JOB_NAME2] [...]

Options:
  -h, --help         print this help and exit

"""

import pickle
import sys
from os import path

# faps is in the parent directory of this script
sys.path.append(path.dirname(path.dirname(path.realpath(__file__))))
# Need the PyNiss and the rest as they are being unpickled
from faps import PyNiss, Structure, Symmetry, Cell, Atom, Guest
from function_switch import to_cif


def niss_to_cif(job_name):
    """Load a simulation Niss and dump out the cif of the current structure."""
    # We just need the job name so don't overcomplicate with configparser
    if job_name[-5:] == ".niss":
        job_name = job_name[:-5]

    # put the niss part back on, even if we took it off already
    niss_name = "%s.niss" % job_name

    if path.exists(niss_name):
        print("Found niss file: %s" % niss_name)
        load_niss = open(niss_name, 'rb')
        my_simulation = pickle.load(load_niss)
        load_niss.close()
    else:
        print("ERROR! %s not found" % niss_name)
        return

    # the method from fapswitch accepts separate components
    cryst = my_simulation.structure
    cif_lines = to_cif(cryst.atoms, cryst.cell, cryst.bonds, job_name)

    # Can't change the name for now; overwrites existing
    cif_file_name = "%s.out.cif" % job_name

    # Write the cif to the file
    output_file = open(cif_file_name, 'w')
    output_file.writelines(cif_lines)
    output_file.close()
    print("Wrote %s" % cif_file_name)


def main():
    """Process multiple niss files given on the command line."""
    # Rudimentary help option
    if "-h" in sys.argv or "--help" in sys.argv:
        print(__doc__)
        raise SystemExit
    # Otherwise everythin is just a job name
    for job_name in sys.argv[1:]:
        niss_to_cif(job_name)


if __name__ == '__main__':

    main()
