"""
configuration for PyTurds

Provides the Options class that will set up directories, defaults, user options
and job options

"""

__all__ = ['Options']

import sys
import os
import ConfigParser
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
        self.args = []
        self.options = {}
        self.load_site_defaults()
        self.load_job_defaults()
        self.commandline()

    def __getattr__(self, item):
        """Maps values to attributes from sources with different priorities."""
        try:
            return self.options.__getitem__(item)
        except KeyError:
            raise AttributeError(item)

    def commandline(self):
        """Specified options, highest priority."""
        usage = "usage: %prog [options] COMMAND"
        parser = OptionParser(usage=usage, version="%prog 0.1",
                              description=__doc__)
        parser.add_option("-v", "--verbose", action="store_true",
                          dest="verbose", default=True,
                          help="write extra information to stdout [default]")
        parser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                          help="silence all output")
        (local_options, local_args) = parser.parse_args()

        #if len(local_args) != 1:
        #    parser.error("No command given (try %prog --help)")

        self.options = local_options
        self.args = local_args

    def load_site_defaults(self):
        """Find where the script is and load defaults"""
        if __name__ != '__main__':
            self.script_dir = os.path.dirname(__file__)
        else:
            self.script_dir = os.path.abspath(sys.path[0])
        site_ini_path = os.path.join(self.script_dir, 'site.ini')
        self.site_ini = ConfigParser.SafeConfigParser()
        self.site_ini.read(site_ini_path)

    def load_job_defaults(self):
        """Find where the job is running and load defaults"""
        self.cwd = os.getcwd()
        job_ini_path = os.path.join(self.cwd, 'job.ini')
        self.job_ini = ConfigParser.SafeConfigParser()
        self.job_ini.read(job_ini_path)
