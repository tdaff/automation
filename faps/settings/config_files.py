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
# Python 3 fix
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from collections import namedtuple
from logging import info, debug, error

DefaultOption = namedtuple('DefaultOption', ['key', 'value', 'type', 'help'])

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
                config_files.append(os.path.join(dot_faps, job_fap))

        all_configs = DictConfigParser()
        # Returns a list of everything that was found
        parsed_configs = all_configs.read(config_files)
        info("Combined settings files".format(", ".join(parsed_configs)))
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


    def get(self, item):
        """Map values from different sources based on priorities."""
        if item in self.__dict__:
            # Instance attributes, such as job_name and job_dir
            debug("an attribute: %s" % item)
            return object.__getattribute__(self, item)
        elif self.options.__dict__.get(item) is not None:
            # Commandline options from optparse where option is set
            debug("an option: %s" % item)
            return self.options.__dict__[item]
        elif item in self.cmdopts:
            # Commandline -o custom key=value options
            debug("a custom -o option: %s" % item)
            return self.cmdopts[item]
        elif self.job_ini.has_option('job_config', item):
            # jobname.fap per-job setings
            debug("a job option: %s" % item)
            return self.job_ini.get('job_config', item)
        elif self.job_type_ini.has_option('job_type', item):
            debug("a job_type option: %s" % item)
            return self.job_type_ini.get('job_type', item)
        elif self.site_ini.has_option('site_config', item):
            debug("a site option: %s" % item)
            return self.site_ini.get('site_config', item)
        elif self.defaults.has_option('defaults', item):
            debug("a default: %s" % item)
            return self.defaults.get('defaults', item)
        else:
            # Most things have a default, but not always. Error properly.
            debug("unspecified option: %s" % item)
            raise AttributeError(item)

    def getbool(self, item):
        """
        Parse option and if the value of item is not already a bool return
        True for "1", "yes", "true" and "on" and False for "0", "no", "false"
        and "off". Case-insensitive.

        """
        value = self.get(item)
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

    def getint(self, item):
        """Return item's value as an integer."""
        value = self.get(item)
        return int(value)

    def getfloat(self, item):
        """Return item's value as a float."""
        value = self.get(item)
        return float(value)

    def gettuple(self, item, dtype=None):
        """Return item's value interpreted as a tuple of 'dtype' [strings]."""
        value = self.get(item)
        # Regex strips bracketing so can't nest, but safer than eval
        value = [x for x in re.split('[\s,\(\)\[\]]*', value) if x]
        if dtype is not None:
            return tuple([dtype(x) for x in value])
        else:
            return tuple(value)


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
