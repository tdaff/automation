#!/usr/bin/env python


import re
import sys
from os import path


class Option(object):
    """Structure to hold information on an individual option."""
    def __init__(self):
        """Put defaults for things like advanced."""
        self.option_name = ""
        self.default = ""
        self.help_text = ""
        self.option_type = ""
        self.advanced = False
        self.section = ""


def parse_defaults():
    """Read all the options from the defaults.ini. Return as a dict."""
    defaults = open(path.join('..', 'defaults.ini'))

    for _header in range(12):
        defaults.readline()

    options = []
    current_option = Option()

    for line in defaults:
        line = line.strip()
        if not line:
            # blank line after, ignore
            continue
        elif line[0] in ['#', ';']:
            line = line.lstrip('#; \t')
            if line.startswith('@type'):
                current_option.option_type = line.split(None, 1)[1]
            elif line.startswith('@section'):
                current_option.section = line.split(None, 1)[1]
            elif line.startswith('@advanced'):
                current_option.advanced = True
            else:
                current_option.help_text += "  %s\n" % line.strip()
        else:
            # Names are simple
            option_name, default = line.split('=')
            current_option.option_name = option_name.strip()

            # Deal with specific types, bools can be real bools
            if current_option.option_type == 'bool':
                if default.strip().lower() in ["1", "yes", "true", "on"]:
                    current_option.default = True
                elif default.strip().lower() in ["0", "no", "false", "off"]:
                    current_option.default = False
            elif 'enum' in current_option.option_type:
                enum_values = re.split(r'[\s,{}]*', current_option.option_type)
                enum_values = [x for x in enum_values if x]
                current_option.option_type = tuple(enum_values)
                current_option.default = default.strip()
            else:
                current_option.default = default.strip()

            options.append(current_option)

            # start fresh for next option
            current_option = Option()

    return options


def main():
    """Parse defaults and generate documentation."""
    rst = []

    for option in parse_defaults():
        rst.append("\n.. envvar:: %s\n\n" % option.option_name)
        # Sphinx doesn't like lots of ````````
        default = "``%s``" % option.default if option.default else ""
        rst.append("  Default: %s\n\n" % default)
        # Appends section description to the option description
        advanced = ":advanced" if option.advanced else ""
        section = "%s%s" % (option.section, advanced)
        rst.append("  **[%s]**\n%s\n" % (section, option.help_text))

    with open('options.inc', 'w') as options_inc:
        options_inc.writelines(rst)


if __name__ == '__main__':
    main()