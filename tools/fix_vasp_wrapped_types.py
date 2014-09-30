#!/usr/bin/env python

"""
fix_vasp_wrapped_types.py

Usage: fix_vasp_wrapped_types.py [filename] [filename] [...]

Faps can create POSCAR files with more atom types on a line than the
outputs can handle and they end up getting wrapped. Use this script
to put the types back on the line so it can be read by other programs.

If filenames are given, only those files will be fixed, otherwise the
CONTCAR CHGCAR LOCPOT files are fixed.

Options:
    -h, --help         display help and exit

"""


import sys

from faps import fix_vasp_wrapped_types


def main():
    """
    Run the fixer. Read filenames from argv and run them through the function.

    """

    if '-h' in sys.argv or '--help' in sys.argv:
        print(__doc__)
        raise SystemExit

    if len(sys.argv) > 1:
        filenames = sys.argv[1:]
    else:
        filenames = ["CONTCAR", "CHGCAR", "LOCPOT"]

    for filename in filenames:
        try:
            if fix_vasp_wrapped_types(filename):
                # Print message only when files have been fixed
                print("Fixed file {}".format(filename))
        except IOError:
            print("Error reading file {}".format(filename))
        except IndexError:
            print("Couldn't convert file {}".format(filename))


if __name__ == '__main__':

    main()
