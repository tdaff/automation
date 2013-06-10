#!usr/bin/env python

"""
Testing faps.util.text functionality

"""

import sys
import unittest

# This is where faps gets imported
sys.path.append('..')
from faps.util.text import unique, ufloat, gfloat, try_int, strip_blanks

class TestUFloat(unittest.TestCase):
    """Uncertainty in floats is ignored."""
    def test_with_uncertainty(self):
        """Typical cases from cif files."""
        self.assertEqual(ufloat('1.03(209)'), 1.03)
        self.assertEqual(ufloat('0(0)'), 0)
        self.assertEqual(ufloat('1982.453(0)'), 1982.453)
        self.assertEqual(ufloat('0(34.2)'), 0)

    def test_no_uncertainty(self):
        """Plain floats with no brackets."""
        self.assertEqual(ufloat('1.03'), 1.03)
        self.assertEqual(ufloat('0'), 0)
        self.assertEqual(ufloat('1982.453'), 1982.453)

    def test_bad_cases(self):
        """Throw error if something bad is given."""
        self.assertRaises(ValueError, ufloat, 'spoon')

class TestGFloat(unittest.TestCase):
    """Gulp outputs fractions somethimes."""
    def test_with_fraction(self):
        """Typical cases from gout files."""
        self.assertEqual(gfloat('1/3'), 1.0/3.0)
        self.assertEqual(gfloat('0/4'), 0)
        self.assertEqual(gfloat('12/11'), 12.0/11.0)

    def test_no_fraction(self):
        """Plain floats with no fraction."""
        self.assertEqual(gfloat('1.03'), 1.03)
        self.assertEqual(gfloat('0'), 0)
        self.assertEqual(gfloat('1982.453'), 1982.453)

    def test_bad_cases(self):
        self.assertRaises(ZeroDivisionError, gfloat, '1/0')
        self.assertRaises(ValueError, gfloat, '1982.453(0)')

class TestTryInt(unittest.TestCase):

    def test_int_values(self):
        """Typical cases from gout files."""
        self.assertEqual(try_int('1'), 1)
        self.assertEqual(try_int('0'), 0)
        self.assertEqual(try_int('34'), 34)
        self.assertEqual(try_int(8.4), 8)

    def test_not_int(self):
        """Plain floats with no fraction."""
        self.assertEqual(try_int('1.03'), 0)
        self.assertEqual(try_int('spoon'), 0)
        self.assertEqual(try_int(''), 0)

    def test_default_value(self):
        self.assertEqual(try_int('1.03', 5), 5)
        self.assertEqual(try_int('spoon', default=7), 7)

if __name__ == '__main__':
    unittest.main(verbosity=3)
