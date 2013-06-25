#!/usr/bin/env python

"""
Symmetry object that stores a symmetry, can operate on positions to generate a
filled cell.

"""

import re


class Symmetry(object):
    """Apply symmetry operations to atomic coordinates."""
    def __init__(self, blob):
        """Read the operation from the argument."""
        self.sym_ops = []
        self.blob = blob
        self.parse_blob()
        self.cell = None

    def parse_blob(self):
        """Interpret a symmetry line from a cif."""
        # convert integers to floats to avoid integer division
        self.sym_ops = [re.sub(r'([\d]+)', r'\1.0', x.strip())
                        for x in re.split(',', self.blob) if x.strip()]

    def trans_frac(self, pos):
        """Apply symmetry operation to the supplied position."""
        new_pos = [eval(sym_op.replace('x', str(pos[0]))
                        .replace('y', str(pos[1]))
                        .replace('z', str(pos[2]))) for sym_op in self.sym_ops]
        # TODO(tdaff): should we translate into cell?
        new_pos = [x%1.0 for x in new_pos]
        return new_pos

