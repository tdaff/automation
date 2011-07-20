"""
configuration for PyFaps

Provides the Options class that will set up directories, defaults, user options
and job options.

"""

__all__ = ['Options']

import sys
import os
import ConfigParser
import __main__
from optparse import OptionParser


class Options(object):
    """
    Transparent options handling.

    A single unified way of dealing with input files and command line options
    delivering sensible defaults for unspecified values. Access options with
    the .get() method (or the method that defines the type). It is recommended
    to create a new instance each time the script is run, otherwise commandline
    options or changed input files will not be picked up.

    """
    def __init__(self):
        """Initialize options from all .ini files and commandline."""
        # use .get*() to read attributes, only access args directly
        self.cwd = ''
        self.script_dir = ''
        self.job_name = ''
        self.args = []
        self.options = {}
        self.cmdopts = {}
        self.defaults = ConfigParser.SafeConfigParser()
        self.site_ini = ConfigParser.SafeConfigParser()
        self.job_ini = ConfigParser.SafeConfigParser()
        # populate options
        self._init_paths()
        self.load_defaults()
        self.load_site_defaults()
        self.load_job_defaults()
        self.commandline()

    def get(self, item):
        """Map values from different sources, based on priorities."""
        # The original
        if item in self.__dict__:
            # Instance attributes, such as cwd
            debug("an attribute: %s" % item)
            return object.__getattribute__(self, item)
        elif self.options.__dict__.get(item) is not None:
            # Commandline options from optparse
            debug("an option: %s" % item)
            return self.options.__dict__[item]
        elif item in self.cmdopts:
            # Commandline -o custom options
            debug("a custom -o option: %s" % item)
            return self.cmdopts[item]
        elif self.job_ini.has_option('job_config', item):
            debug("a job option: %s" % item)
            return self.job_ini.get('job_config', item)
        elif self.site_ini.has_option('site_config', item):
            debug("a site option: %s" % item)
            return self.site_ini.get('site_config', item)
        elif self.defaults.has_option('defaults', item):
            debug("a default: %s" % item)
            return self.defaults.get('defaults', item)
        else:
            # Shouldn't get here; everything should have a default!
            print "unspecified option: %s" % item
            raise AttributeError(item)

    def getbool(self, item):
        """
        Parse option and if item is not already a bool return True for "1",
        "yes", "true" and "on" and False for "0", "no", "false", and "off".
        Case-insensitive.

        """
        value = self.get(item)
        if isinstance(value, bool):
            return value
        elif isinstance(value, basestring):
            if value.lower() in ["1", "yes", "true", "on"]:
                return True
            elif value.lower() in ["0", "no", "false", "off"]:
                return False
            else:
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

    def gettuple(self, item):
        """Return item's value interpreted as a tuple (best guess)."""
        value = self.get(item)
        if isinstance(item, basestring):
            # NOTE: be careful, eval can be dangerous
            # TODO(tdaff): catch invalid values
            value = eval(value, {}, {})
            return tuple(value)
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
        self.cwd = os.getcwd()


    def commandline(self):
        """Specified options, highest priority."""
        usage = "usage: %prog [options] [COMMAND] JOB_NAME"
        # use description for the script, not for this module
        parser = OptionParser(usage=usage, version="%prog 0.1",
                              description=__main__.__doc__)
        parser.add_option("-v", "--verbose", action="store_true",
                          dest="verbose", default=True,
                          help="write extra information to stdout [default]")
        parser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                          help="silence all output")
        parser.add_option("-o", "--option", action="append", dest="cmdopts",
                          help="set custom options as key=value pairs")
        parser.add_option("-i", "--interactive", action="store_true",
                          dest="interactive", help="enter interactive mode")
        parser.add_option("-m", "--import", action="store_true",
                          dest="import", help="try and import old data")
        parser.add_option("-n", "--nosub", action="store_true", dest="nosub",
                          help="create input files only, do not run any jobs")
        (local_options, local_args) = parser.parse_args()

        if len(local_args) == 0:
            parser.error("No arguments given (try %prog --help)")
        elif len(local_args) == 1:
            # continue a full job if just the job name is given
            self.job_name = local_args.pop()
            local_args = ['run']
        else:
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
        """Load program defaults"""
        default_ini_path = os.path.join(self.script_dir, 'defaults.ini')
        self.defaults.read(default_ini_path)

    def load_site_defaults(self):
        """Find where the script is and load defaults"""
        site_ini_path = os.path.join(self.script_dir, 'site.ini')
        self.site_ini.read(site_ini_path)

    def load_job_defaults(self):
        """Find where the job is running and load defaults"""
        job_ini_path = os.path.join(self.cwd, 'job.ini')
        self.job_ini.read(job_ini_path)


def debug(msg):
    """Print for debugging statements."""
    print("DEBUG: %s" % msg)



def options_test():
    """Try and read a few options from different sources."""
    testopts = Options()
    print(testopts.get('job_name'))
    print(testopts.get('cmdopts'))
    print(testopts.get('args'))
    print(testopts.get('verbose'))
    print(testopts.get('script_dir'))
    print(testopts.get('interactive'))
    print(testopts.get('whot'))
    print(testopts.get('repeat_exe'))


if __name__ == '__main__':
    options_test()
