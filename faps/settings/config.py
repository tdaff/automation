#!/usr/bin/env python

"""
configuration for faps

Provides the Options class that will transparently handle the different option
sources through the .get() method. Pulls in defaults, site and job options plus
command line customisation. Instantiating Options will set up the logging for
the particular job.

"""

__all__ = ['Options']

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

import __main__

job_path = os.getcwd()

class Options(object):
    """
    Transparent options handling.

    A single unified way of dealing with input files and command line options
    delivering sensible defaults for unspecified values. Access options with
    the .get() method, or the method that specifies the expected type. It is
    recommended to replace with a new instance each time the script is run,
    otherwise commandline options or changed input files will not be picked up.

    """
    def __init__(self, job_name=None):
        """Initialize options from all .ini files and the commandline."""
        # use .get{type}() to read attributes, only access args directly
        self.job_dir = ''
        self.script_dir = ''
        self.job_name = job_name
        self.args = []
        self.options = {}
        self.cmdopts = {}
        self.defaults = configparser.SafeConfigParser()
        self.site_ini = configparser.SafeConfigParser()
        self.job_ini = configparser.SafeConfigParser()
        # populate options
        self._init_paths()
        self.commandline()
        self._init_logging()
        self.load_defaults()
        self.load_site_defaults()
        self.load_job_defaults()
        if self.options.job_type:
            self.job_type_ini = configparser.SafeConfigParser()
            self.load_job_type(self.options.job_type)
        else:
            self.job_type_ini = NullConfigParser()

    def get(self, item):
        """Map values from different sources based on priorities."""
        if item in self.__dict__:
            # Instance attributes, such as job_name and job_dir
            debug("an attribute: %s" % item)
            return object.__getattribute__(self, item)
        elif self.options.__dict__.get(item) is not None:
            # Commandline options from optparse where option is set
            debug("an option: %s" % item)
            return self.options.__dict__[item]
        elif item in self.cmdopts:
            # Commandline -o custom key=value options
            debug("a custom -o option: %s" % item)
            return self.cmdopts[item]
        elif self.job_ini.has_option('job_config', item):
            # jobname.fap per-job setings
            debug("a job option: %s" % item)
            return self.job_ini.get('job_config', item)
        elif self.job_type_ini.has_option('job_type', item):
            debug("a job_type option: %s" % item)
            return self.job_type_ini.get('job_type', item)
        elif self.site_ini.has_option('site_config', item):
            debug("a site option: %s" % item)
            return self.site_ini.get('site_config', item)
        elif self.defaults.has_option('defaults', item):
            debug("a default: %s" % item)
            return self.defaults.get('defaults', item)
        else:
            # Most things have a default, but not always. Error properly.
            debug("unspecified option: %s" % item)
            raise AttributeError(item)

    def getbool(self, item):
        """
        Parse option and if the value of item is not already a bool return
        True for "1", "yes", "true" and "on" and False for "0", "no", "false"
        and "off". Case-insensitive.

        """
        value = self.get(item)
        if isinstance(value, bool):
            return value
        # Can't use isinstance with basestring to be 2.x and 3.x compatible
        # fudge it by assuming strings can be lowered
        elif hasattr(value, 'lower'):
            if value.lower() in ["1", "yes", "true", "on"]:
                return True
            elif value.lower() in ["0", "no", "false", "off"]:
                return False
            else:
                # Not a valid bool
                raise ValueError(value)
        else:
            return bool(item)

    def getint(self, item):
        """Return item's value as an integer."""
        value = self.get(item)
        return int(value)

    def getfloat(self, item):
        """Return item's value as a float."""
        value = self.get(item)
        return float(value)

    def gettuple(self, item, dtype=None):
        """Return item's value interpreted as a tuple of 'dtype' [strings]."""
        value = self.get(item)
        # Regex strips bracketing so can't nest, but safer than eval
        value = [x for x in re.split('[\s,\(\)\[\]]*', value) if x]
        if dtype is not None:
            return tuple([dtype(x) for x in value])
        else:
            return tuple(value)

    def _init_paths(self):
        """Find the script directory and set up working directory"""
        # Where the script is has the config defaults.
        if __name__ != '__main__':
            self.script_dir = os.path.dirname(__file__)
        else:
            self.script_dir = os.path.abspath(sys.path[0])
        # Where we run the job.
        self.job_dir = os.getcwd()


    def commandline(self):
        """Specified options, highest priority."""
        usage = "usage: %prog [options] [COMMAND] JOB_NAME"
        # use description for the script, not for this module
        parser = OptionParser(usage=usage, version="%prog 0.1",
                              description=__main__.__doc__)
        parser.add_option("-v", "--verbose", action="store_true",
                          dest="verbose",
                          help="output extra debugging information")
        parser.add_option("-q", "--quiet", action="store_true",
                          dest="quiet", help="only output warnings and errors")
        parser.add_option("-s", "--silent", action="store_true",
                          dest="silent", help="no terminal output")
        parser.add_option("-p", "--plain", action="store_true",
                          dest="plain", help="do not colourise or wrap output")
        parser.add_option("-o", "--option", action="append", dest="cmdopts",
                          help="set custom options as key=value pairs")
        parser.add_option("-i", "--interactive", action="store_true",
                          dest="interactive", help="enter interactive mode")
        parser.add_option("-m", "--import", action="store_true",
                          dest="import", help="try and import old data")
        parser.add_option("-n", "--no-submit", action="store_true",
                          dest="no_submit",
                          help="create input files only, do not run any jobs")
        parser.add_option("-j", "--job-type", dest="job_type",
                          help="user preconfigured job settings")
        parser.add_option("-d", "--daemon", action="store_true", dest="daemon",
                          help="run [lube] as a server and await input")
        (local_options, local_args) = parser.parse_args()

        # job_name may or may not be passed or set initially
        if self.job_name:
            if self.job_name in local_args:
                local_args.remove(self.job_name)
        elif len(local_args) == 0:
            parser.error("No arguments given (try %prog --help)")
        else:
            # Take the last argument as the job name
            self.job_name = local_args.pop()

        # key value options from the command line
        if local_options.cmdopts is not None:
            for pair in local_options.cmdopts:
                if '=' in pair:
                    pair = pair.split('=')
                    self.cmdopts[pair[0]] = pair[1]
                else:
                    self.cmdopts[pair] = True

        self.options = local_options
        # Args are only the COMMANDS for the run
        self.args = [arg.lower() for arg in local_args]

    def load_defaults(self):
        """Load program defaults."""
        # ConfigParser requires header sections so we add them to a StringIO
        # of the file if they are missing. 2to3 should also deal with the
        # renamed modules.
        default_ini_path = os.path.join(self.script_dir, 'defaults.ini')
        try:
            filetemp = open(default_ini_path, 'r')
            default_ini = filetemp.read()
            filetemp.close()
            if not '[defaults]' in default_ini.lower():
                default_ini = '[defaults]\n' + default_ini
            default_ini = StringIO(default_ini)
        except IOError:
            # file does not exist so we just use a blank string
            debug('Default options not found! Something is very wrong.')
            default_ini = StringIO('[defaults]\n')
        self.defaults.readfp(default_ini)

    def load_site_defaults(self):
        """Find where the script is and load defaults"""
        site_ini_path = os.path.join(self.script_dir, 'site.ini')
        try:
            filetemp = open(site_ini_path, 'r')
            site_ini = filetemp.read()
            filetemp.close()
            if not '[site_config]' in site_ini.lower():
                site_ini = '[site_config]\n' + site_ini
            site_ini = StringIO(site_ini)
        except IOError:
            # file does not exist so we just use a blank string
            debug("No site options found; using defaults")
            site_ini = StringIO('[site_config]\n')
        self.site_ini.readfp(site_ini)

    def load_job_defaults(self):
        """Find where the job is running and load defaults"""
        job_ini_path = os.path.join(self.job_dir, self.job_name + '.fap')
        try:
            filetemp = open(job_ini_path, 'r')
            job_ini = filetemp.read()
            filetemp.close()
            if not '[job_config]' in job_ini.lower():
                job_ini = '[job_config]\n' + job_ini
            job_ini = StringIO(job_ini)
            debug("Job options read from %s" % job_ini_path)
        except IOError:
            # file does not exist so we just use a blank string
            debug("No job options found; using defaults")
            job_ini = StringIO('[job_config]\n')
        self.job_ini.readfp(job_ini)

    def load_job_type(self, job_type):
        """Find where the job is running and load defaults"""
        home_dir = os.path.expanduser('~')
        job_type_ini_path = os.path.join(home_dir, '.faps', job_type + '.fap')
        try:
            filetemp = open(job_type_ini_path, 'r')
            job_type_ini = filetemp.read()
            filetemp.close()
            if not '[job_type]' in job_type_ini.lower():
                job_type_ini = '[job_type]\n' + job_type_ini
            job_type_ini = StringIO(job_type_ini)
            debug("Job type options read from %s" % job_type_ini_path)
        except IOError:
            # file does not exist so we just use a blank string
            error("Job type '%s' specified but options file '%s' not found" %
                  (job_type, job_type_ini_path))
            job_type_ini = StringIO('[job_config]\n')
        self.job_type_ini.readfp(job_type_ini)



class NullConfigParser(object):
    """Use in place of a blank ConfigParser that has no options."""
    def __init__(self, *args, **kwargs):
        """This is empty, so do nothing."""
        pass

    def has_option(*args, **kwargs):
        """Always return Fasle as there are no options."""
        return False


if __name__ == '__main__':
    options_test()
