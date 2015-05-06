#!/usr/bin/env python

import sys

defaults = open('defaults.ini')

for _header in xrange(12):
    defaults.readline()

rst = []

this_option = []

for line in defaults:
    line = line.strip().replace('*', '\\*')
    if not line:
        continue
    elif line[0] in ['#', ';']:
        this_option.append("  %s\n" % line.lstrip('#; \t'))
    else:
        option, default = line.split('=')
        rst.append("\n.. envvar:: %s\n\n" % option.strip())
        rst.append("  Default: %s\n\n" % default.strip())
        rst.extend(this_option)
        this_option = []

sys.stdout.writelines(rst)
