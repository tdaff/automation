#!/usr/bin/env python

"""
Small system and filesystem processing utilities.

"""

import glob
import os
import shutil
import subprocess
import sys
from logging import warning, debug, error, info, critical
from os import path


def mkdirs(directory):
    """Create a directory if it does not exist."""
    if not path.exists(directory):
        os.makedirs(directory)


def terminate(exit_code=0):
    """Exit and announce if faps is terminating normally (default)."""
    if exit_code == 0:
        info("Faps terminated normally")
        raise SystemExit
    else:
        warning("Abnormal termination of faps; exit code %i" % exit_code)
        raise SystemExit(exit_code)


def move_and_overwrite(src, dest):
    """Move src to dest and overwrite if it is an existing file."""
    # As src and dest can be files or directories, do some checks.
    if path.exists(dest):
        if path.isdir(dest):
            dest_full = path.join(dest, path.basename(src))
            if path.exists(dest_full):
                if path.isfile(dest_full):
                    os.remove(dest_full)
                    shutil.move(src, dest)
                else:
                    raise OSError("Directory %s already exists" % dest_full)
            else:
                shutil.move(src, dest)
        elif path.isfile(dest):
            os.remove(dest)
            shutil.move(src, dest)
        else:
            raise OSError("%s is not a folder or file" % dest)
    else:
        shutil.move(src, dest)


def try_symlink(src, dest):
    """Delete an existing dest file, symlink a new one if possible."""
    if path.lexists(dest):
        os.remove(dest)
    try:
        os.symlink(src, dest)
    except AttributeError:
        shutil.copy(src, dest)


def sys_argv_strip(argument):
    """Remove an argument from the sys.argv if it is there."""
    while argument in sys.argv:
        sys.argv.remove(argument)
    return 0


def remove_files(files, directory='.'):
    """
    Delete any of the files if they exist using standard globbing,
    remove empty directories, or ignore silently if not found.

    """
    del_list = []
    for file_name in files:
        del_list.extend(glob.glob(path.join(directory, file_name)))
    for del_name in del_list:
        debug("deleting %s" % del_name)
        try:
            os.remove(del_name)
        except OSError:
            try:
                os.rmdir(del_name)
            except OSError:
                debug("Directory not empty: %s?" % del_name)


def compress_files(files, directory='.'):
    """Gzip any big files to keep."""
    zip_list = []
    for file_name in files:
        zip_list.extend(glob.glob(path.join(directory, file_name)))
    for zip_name in zip_list:
        debug("compressing %s" % zip_name)
        gzip_command = ['gzip', '-f', zip_name]
        subprocess.call(gzip_command)
        #TODO(tdaff): internal gzip

