#!/usr/bin/env python

"""
PyTurds

Automated high throughput strucutre analysis
"""

from optparse import OptionParser
import pickle


class Options(object):
    """
    A single object to deal with input files and command line options but
    delivering reasonable defaults for unspecified values.
    """
    def __init__(self):
        self._commandline()

    def _commandline(self):
        """Define usage options, some error checking"""
        usage = "usage: %prog [options] COMMAND"
        parser = OptionParser(usage=usage, version="%prog 0.1",
                              description=__doc__)
        parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                          help="write extra information to stdout [default]",
                          default=True)
        parser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                          help="silence all output")
        (local_options, local_args) = parser.parse_args()

        #if len(local_args) != 1:
        #    parser.error("No command given (try %prog --help)")

        self.options = local_options
        self.args = local_args

class Simulation(object):
    """
    A single set of calculations for a structure
    """
    def __init__(self, options):
        self.options = options




if __name__ == '__main__':
    global_options = Options()
    my_simulation = Simulation(global_options)
    print(pickle.dumps(my_simulation))
