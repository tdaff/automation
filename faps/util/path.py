"""
Job handler

Machine specific job submission and tracking routines. Provides a unified
interface for all job handling. Implements the JobHandler class which will be
initialized to the machine the calculations are running on from the provided
options.

"""

import os
import re
import shlex
import subprocess
import sys
import time
from logging import debug, info, warn, error
from os import path
from subprocess import Popen, PIPE, STDOUT

def argstrip(arglist):
    """Some options might be best removed before resubmission."""
    to_remove = ['-i', '--interactive', '-m', '--import']
    newargs = list(arglist)
    for item in to_remove:
        while item in newargs:
            newargs.remove(item)
    return newargs


def check_program(program, mpi=False):
    """Test to see if the exe is in the path and might work."""
    exe = which(program)
    if exe is None:
        error("Could not find %s in path, job will probably fail." % program)
        return False
    elif mpi:
        try:
            binary = open(exe, 'rb')
            if re.search('mpi_init', binary.read(), re.IGNORECASE):
                return True
            else:
                warn("%s doesn't appear to be an mpi executable." % program)
                return False
        except IOError:
            return False
    else:
        # TODO(tdaff): No easy way to check for OMP in general?
        return True


def is_exe(fpath):
    """Return executability of a file."""
    return path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """Return the equivalent of the 'which' command."""

    fpath, _fname = path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for env_path in os.environ["PATH"].split(os.pathsep):
            exe_file = path.join(env_path, program)
            if is_exe(exe_file):
                return exe_file

    return None
