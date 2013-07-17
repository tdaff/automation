#!/usr/bin/env python

"""
Output logging setup for faps. If this is not imported before logging
begins weird things will happen.

"""

import copy
import logging
import sys


def init_logging(basename, verbosity=0):
    """
    Setup the logging to terminal and .flog file, with levels as required.
    Packages are expected to import the provided logger and functions rather
    than roll their own. The verbose count should be how many times -v or
    -q appear on the command line.

    """


    root_logger = logging.getLogger()

    # Have the logger itself set with the lowest possible level
    root_logger.setLevel(logging.DEBUG)
    # Reset any handlers that might have been set accidentally
    root_logger.handlers = []

    # Always at least INFO in .flog
    file_level = logging.INFO
    flog_filename = '{}.flog'.format(basename)

    if verbosity <= -2:
        # -qq
        stdout_level = logging.CRITICAL
    elif verbosity <= -1:
        # -q
        stdout_level = logging.ERROR
    elif verbosity >= 1:
        # -v
        stdout_level = logging.DEBUG
        file_level = logging.DEBUG
    else:
        stdout_level = logging.INFO

    # Easier to do simple file configuration then add the stdout
    file_handler = logging.FileHandler(flog_filename)
    file_handler.setLevel(file_level)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s %(message)s',
                                  datefmt='%Y%m%d %H:%M:%S')
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

    # Make these uniform widths
    logging.addLevelName(10, '--')
    logging.addLevelName(20, '>>')
    logging.addLevelName(30, '**')
    logging.addLevelName(40, '!!')
    logging.addLevelName(50, 'XX')

    # Use nice coloured console output
    console = ColouredConsoleHandler(stream=sys.stdout)
    console.setLevel(stdout_level)
    formatter = logging.Formatter('%(levelname)s %(message)s')
    console.setFormatter(formatter)
    # add the handler to the root logger
    root_logger.addHandler(console)


class ColouredConsoleHandler(logging.StreamHandler):
    """Makes colourised output for the console."""
    def emit(self, record):
        """Colourise leve id and emit a record."""
        # Need to make a actual copy of the record
        # to prevent altering the message for other loggers
        myrecord = copy.copy(record)
        levelno = myrecord.levelno
        if levelno >= 50:  # CRITICAL / FATAL
            front = '\033[30;41m'  # black/red
        elif levelno >= 40:  # ERROR
            front = '\033[30;41m'  # black/red
        elif levelno >= 30:  # WARNING
            front = '\033[30;43m'  # black/yellow
        elif levelno >= 20:  # INFO
            front = '\033[30;42m'  # black/green
        elif levelno >= 10:  # DEBUG
            front = '\033[30;46m'  # black/cyan
        else:  # NOTSET and anything else
            front = '\033[0m'  # normal

        myrecord.levelname = '{}{}\033[0m'.format(front, myrecord.levelname)
        logging.StreamHandler.emit(self, myrecord)
