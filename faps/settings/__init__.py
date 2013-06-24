"""
Faps settings and configurations.

"""

import os
from os import path


# Paths that define where to find things in the code
# these are not easy to work out within the code so
# consolidate them here.
settings_path = os.path.dirname(os.path.abspath(__file__))
data_path = path.join(settings_path, '..', 'data')
dot_faps_path = path.join(path.expanduser('~'), '.faps')