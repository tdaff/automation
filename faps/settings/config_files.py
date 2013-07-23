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
    from configparser import SafeConfigParser
except ImportError:
    from ConfigParser import SafeConfigParser
import copy
import logging
import os
import re
import sys
import textwrap

from collections import namedtuple
from logging import info, debug, error

DefaultOption = namedtuple('DefaultOption', ['key', 'value', 'type', 'help'])

# Use in the case of unspecified defaults
_default = DefaultOption(key='', value='', type='str', help='Empty option')


class Options(object):
    """
    Transparent options handling.

    A single unified way of dealing with input files and command line options,
    delivering sensible defaults for unspecified values. Access options with
    the .get() method which are typed by the default setting to the
    expected type.

    All components must interface to .add_defaults() to register their
    expected types and fallbacks.

    """
    def __init__(self, arguments=None, path=None):
        """
        Initialize options from all config files and the commandline.
        It is expected that arguments has a job_name attribute and path
        contains the paths to look for config files.

        """
        self.job_name = getattr(arguments, 'job_name', 'default')
        self.path = path
        # use .get() to read attributes, only access args directly
        self.defaults = {}
        # get any defaults first

        # combine user options in the order:
        #  -> job types, in order
        #  -> job options
        #  -> command line options

        self.options = {}

        # build a list of all possible config files:
        config_files = []

        # Preset job types pulled in at the bottom
        if hasattr(path, 'settings'):
            presets = os.path.join(path.settings, 'presets')
            for job_type in getattr(arguments, 'job_type', []):
                job_fap = '{}.fap'.format(job_type)
                config_files.append(os.path.join(presets, job_fap))

        # User custom job types
        if hasattr(path, 'dot_faps'):
            dot_faps = path.dot_faps
            for job_type in getattr(arguments, 'job_type', []):
                job_fap = '{}.fap'.format(job_type)
                print(job_fap)
                config_files.append(os.path.join(dot_faps, job_fap))

        all_configs = DictConfigParser()
        # Returns a list of everything that was found
        parsed_configs = all_configs.read(config_files)
        info("Combined settings files {}".format(", ".join(parsed_configs)))
        self.options = all_configs.as_dict()

        # parse and update the options with
        # section.key=value options from the command line
        if hasattr(arguments, 'cmdopts'):
            for cmdopt in arguments.cmdopts:
                try:
                    section, kv_pair = cmdopt.split('.', 1)  # maxsplit=1
                    if '=' in kv_pair:
                        key, value = kv_pair.split('=')
                    else:
                        # Set as true when no value specified
                        key, value = kv_pair, "True"
                    # Section might not exist already
                    if not section in self.options:
                        self.options['section'] = {}
                    # Should be fine to insert now
                    self.options['section'][key] = value
                except ValueError:
                    # Probably no section specified
                    error("Invalid option in command line: {}".format(cmdopt))

        debug("My options look like {}".format(self.options))


    def get(self, section, item):
        """Return the (typed) value for the option."""
        # A lookup for the conversion functions helps deal with strings
        vtypes = {
            'str': str,
            'float': float,
            'int': int,
            'bool': _bool,
            'tuple': _tuple,
        }

        # Make sure that we know the option
        try:
            default = self.defaults[item]
        except KeyError:
            warning('Unknown option {}'.format(item))
            default = _default

        # need to know the expected options type
        vtype = default.type
        if 'tuple' in vtype or 'list' in vtype:
            # split on non word character
            dtype = re.split('\W*', vtype)[1]
            if dtype not in vtypes:
                warning("Unknown tuple dtype {}".format(dtype))
                dtype = 'str'
        elif vtype not in vtypes:
            warning("Unknown option type {}".format(vtype))
            vtype = 'str'

        # Check that options will be okay with the request
        if not section in self.options:
            warning("Unidentified section - shouldn't happen! Using default")
            return vtypes[vtype](default.value)
            #raise AttributeError(item)
        elif item in self.options[section]:
            # Directly specified
            return vtypes[vtype](self.options[section][item])
        elif item in self.options['global']:
            # Specified globally
            return vtypes[vtype](self.options['global'][item])
        else:
            # Just use the default
            return vtypes[vtype](default.value)


    def make_section_getter(self, section):
        """Wrap the getter function for a fixed section of the options."""
        # closured section cleans up the interface a bit
        def section_getter(item):
            """return .get() for the item, expects section from outer scope."""
            return self.get(section, item)


    def register_default(self, key, value, type, help):
        """
        Make a new default available in the options. Required arguments are:

        key:   The option name as it would be in the options file.
        value: The default value of the option, strings
               are fine as they will be converted.
        type:  The datatype expected for the option.
        help:  A description of the option for documentation.

        """
        new_default = DefaultOption(key, value, type, help)
        self.defaults[key] = new_default

    @property
    def targets(self):
        option_targets = set(self.options)
        option_targets.remove('global')
        if len(option_targets) == 0:
            return {'default'}
        else:
            return option_targets


class DictConfigParser(SafeConfigParser):
    """
    Derivative of configparser class that provides dict and some
    extra sectioning.

    """
    def as_dict(self):
        """Return a simple non-ordered dict of all the options."""
        drepr = {}
        for section in self.sections():
            drepr[section] = dict(self.items(section))
        return drepr



def _bool(value):
    """
    Parse the value and return a boolean. Strings (case-insensitive) return
    True for "1", "yes", "true" and "on" and False for "0", "no", "false"
    and "off". Other types return their Python boolean values.

    """
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

def _tuple(value, dtype=None):
    """
    Split value into a tuple, ignoring whitespace and bracketing, and
    optionally apply 'dtype' function to each item [str]."""
    # Regex strips bracketing so can't nest, but safer than eval
    value = [x for x in re.split('[\s,\(\)\[\]]*', value) if x]
    if dtype is not None:
        return tuple([dtype(x) for x in value])
    else:
        return tuple(value)
