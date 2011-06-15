#!/usr/bin/env python

"""
PyTurds

Automated high throughput strucutre analysis
"""

def get_opts():
    """Define usage options, some error checking"""
    usage = "usage: %prog [options] INFILE FX FY FZ OUTPREFIX"
    parser = OptionParser(usage=usage, version="%prog 0.1",
                          description=__doc__)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="write extra information to stdout [default]",
                      default=True)
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                      help="silence all output")
    parser.add_option("-c", "--crop", action="store_true", dest="crop",
                      help="crop atoms to those in the cell [default]",
                      default=True)
    parser.add_option("-n", "--no-crop", action="store_false", dest="crop",
                      help="do not crop atoms; use original configuration")
    parser.add_option("-z", "--zoom", dest="zoom", type='float',
                      help="resize the grid with a spline interpolation")
    parser.add_option("-A", action="store_true", dest="flatten_a",
                      help="flatten along the a axis")
    parser.add_option("-B", action="store_true", dest="flatten_b",
                      help="flatten along the b axis")
    parser.add_option("-C", action="store_true", dest="flatten_c",
                      help="flatten along the c axis")
    (local_options, local_args) = parser.parse_args()

    if len(local_args) != 5:
        parser.error("Five arguments are required (try %prog --help)")

    return (local_options, local_args)
