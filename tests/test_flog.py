#!usr/bin/env python

"""
Testing logging functionality

"""

import os
import logging
import sys
import unittest

# This is where faps gets imported
sys.path.append('..')
from faps.settings.flog import init_logging


class TestFileLogging(unittest.TestCase):

    def setUp(self):
        """Create a verbose logging setup."""
        init_logging('test', 9)

    def test_log_to_file(self):
        """Spew some messages and check they make it to the file."""
        logging.error("Error Message")
        logging.warning("Warning Message")
        logging.info("Info Message")
        logging.debug("Debug Message")
        self.assertTrue(os.path.exists('test.flog'))
        with open('test.flog') as test_flog:
            file_contents = test_flog.read()
        self.assertTrue("Error Message" in file_contents)
        self.assertTrue("Warning Message" in file_contents)
        self.assertTrue("Info Message" in file_contents)
        self.assertTrue("Debug Message" in file_contents)

    def tearDown(self):
        """
        Kill any remaining logging leftovers.
        """
        logging.shutdown()
        try:
            os.remove('test.flog')
        except OSError:
            # File not found
            pass


if __name__ == '__main__':
    unittest.main(verbosity=0, buffer=True)
