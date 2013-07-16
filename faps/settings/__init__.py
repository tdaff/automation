"""
Faps settings and configurations.

The module has a lot of side effects that involve reading configuration
and setting up output streams. It should probably be imported before you
do anything at all.

"""

import os
from os import path

# We can read the commandline and setup logging before everything else

from faps.settings.flog import init_logging
from faps.settings.commandline import commandline as _commandline


commandline_arguments = _commandline()

init_logging(basename=commandline_arguments.job_name,
             verbosity=commandline_arguments.verbosity)


# Paths that define where to find things in the code
# these are not easy to work out within the code so
# consolidate them here.
settings_path = os.path.dirname(os.path.abspath(__file__))
data_path = path.join(settings_path, '..', 'data')
dot_faps_path = path.join(path.expanduser('~'), '.faps')


# Configuration data needs to get collected once everything has been set up
# for the logging and command line

import logging

flog = logging.getLogger('flog')
debug = flog.debug
info = flog.info
warning = flog.warning
error = flog.error
critical = flog.critical

#from faps.settings.config import Options

#config = Options(commandline_arguments)

