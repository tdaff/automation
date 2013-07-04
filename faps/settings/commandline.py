#!/usr/bin/env python

"""
Command line argument parsing.

"""


import argparse

# Python 3 fix
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import copy
import logging
import os
import re
import sys
import textwrap
# Python 3 fix
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from optparse import OptionParser
from logging import debug, error


DESCRIPTION = """
faps -- Frontend for Automated Adsorption Analysis of Porous Solids.

Strucutre adsorption property analysis for high throughput processing.
When run a script, faps will automatically run complete analysis on a 
structure. Sensible defaults are implemented, but calculations can be
easily customised.
"""

# Need to derive a custom argparse action for decrementing values.
# This is needed so that -v and -q can have opposing actions.
class DecreaseAction(argparse.Action):
    """Decrease the destination by 1 for each call."""
    def __call__(self, parser, namespace, values, option_string=None):
        """Decrease destination by one."""
        previous_value = getattr(namespace, self.dest)
        setattr(namespace, self.dest, previous_value-1)


def commandline():
    """Pull everything from the command line."""
    usage = "usage: %prog [options] [COMMAND] JOB_NAME"
    # use description for the script, not for this module
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        dest="verbose", help="Increase verbosity. Specify "
                        "more times for more debugging information. Cancels "
                        "'--quiet'.")
    parser.add_argument("-q", "--quiet", action=DecreaseAction, nargs=0,
                        dest="verbose", help="Decrease verbosity. Specify "
                        "more times for less output. Cancels '--verbose'.")
    parser.add_argument("-o", "--option", action="append", dest="cmdopts",
                        help="Set program options as key=value pairs. Use "
                        "\"quotation marks\" if options contain spaces.")
    parser.add_argument("-j", "--job-type", dest="job_type",
                        help="Read preconfigured job settings from "
                        "job-type.fap in the user ~/.faps/ directory")
    parser.add_argument('job_name', help="Name for job")
    #TODO(tdaff): cleanup or implement    
    #parser.add_option("-i", "--interactive", action="store_true",
    #                    dest="interactive", help="enter interactive mode")
    #parser.add_option("-m", "--import", action="store_true",
    #                    dest="import", help="try and import old data")
    #parser.add_option("-n", "--no-submit", action="store_true",
    #                    dest="no_submit",
    #                    help="create input files only, do not run any jobs")
    #parser.add_option("-d", "--daemon", action="store_true", dest="daemon",
    #                    help="run [lube] as a server and await input")
    local_args = parser.parse_args()

    cmdopts = {}
    # key value options from the command line
    if local_args.cmdopts is not None:
        for pair in local_args.cmdopts:
            if '=' in pair:
                pair = pair.split('=')
                cmdopts[pair[0]] = pair[1]
            else:
                cmdopts[pair] = True
    
    return(local_args, cmdopts)

#    self.options = local_options
    # Args are only the COMMANDS for the run
#    self.args = [arg.lower() for arg in local_args]

print(commandline())