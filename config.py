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

    A single object to deal with input files and command line options but
    delivering reasonable defaults for unspecified values. Accessible through
    attributes of the instance.

    """
    def __init__(self):
        """Initialize options as site.ini, job.ini and commandline."""
        self.cwd = ''
        self.script_dir = ''
        self.job_name = ''
        self.args = []
        self.options = {}
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
        """Map values to attributes from many sources, based on priorities."""
        if item in self.__dict__:
            print "an attribute: %s" % item
            return object.__getattribute__(self, item)
        elif item in self.options.__dict__:
            print "an option: %s" % item
            return self.options.__dict__[item]
        elif self.job_ini.has_option('job_config', item):
            print "a job option: %s" % item
            return self.job_ini.get('job_config', item)
        elif self.site_ini.has_option('site_config', item):
            print "a site option: %s" % item
            return self.site_ini.get('site_config', item)
        elif self.defaults.has_option('defaults', item):
            print "a default: %s" % item
            return self.defaults.get('defaults', item)
        else:
            print "unspecified option: %s" % item
            raise AttributeError(item)

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


if __name__ == '__main__':
    testopts = Options()
    print(testopts.get('job_name'))
    print(testopts.get('cmdopts'))
    print(testopts.get('args'))
    print(testopts.get('verbose'))
    print(testopts.get('script_dir'))
    print(testopts.get('cwdt'))
    print(testopts.get('repeat_exe'))
