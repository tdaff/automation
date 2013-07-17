"""
Faps settings and configurations.

The module has a lot of side effects that involve reading configuration
and setting up output streams. It should probably be imported before you
do anything at all.

"""

import os as _os
from os import path as _path

# Read the commandline first as logging can depend on it

from faps.settings.commandline import commandline as _commandline

_arguments = _commandline()

# Setup logging before everything else so logging functions
# can be used anywhere

from faps.settings.flog import init_logging as _init_logging

_init_logging(basename=_arguments.job_name,
              verbosity=_arguments.verbosity)

import logging as _logging

flog = _logging.getLogger('flog')
debug = flog.debug
info = flog.info
warning = flog.warning
error = flog.error
critical = flog.critical


# Paths that define where to find things in the code
# these are not easy to work out within the code so
# consolidate them here.
class _Path(object):
    pass

path = _Path()
path.start = _os.getcwd()
path.job = _path.join(path.start, '{}_faps'.format(_arguments.job_name))
path.settings = _path.dirname(_path.abspath(__file__))
path.data = _path.join(path.settings, '..', 'data')
path.dot_faps = _path.join(_path.expanduser('~'), '.faps')


# Configuration data needs to get collected once everything has been set up
# for the logging and command line

#from faps.settings.config_files import Options

#config = Options(commandline_arguments)

