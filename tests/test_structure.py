#!usr/bin/env python

"""
Testing Structure functionality

"""

import sys
import unittest

# This is where faps gets imported
sys.path.append('..')
from faps.core.structure import Structure


class TestInitialiseStructure(unittest.TestCase):
    
    def test_initialise_structure(self):
        """Just make a structure."""
        new_structure = Structure()
        self.assertIsInstance(new_structure, Structure)
        self.assertIsNone(new_structure.name)


class TestAttributes(unittest.TestCase):

    def setUp(self):
        """Create empty cell to work with."""
        self.structure = Structure()

    def test_add_and_remove_attributes(self):
        """Add attributes, remove attributes and test that they work"""
        self.assertFalse(hasattr(self.structure, 'test'))
        self.structure.test = 101
        self.assertTrue(hasattr(self.structure, 'test'))
        self.assertEqual(self.structure.test, 101)
        del self.structure.test
        self.assertFalse(hasattr(self.structure, 'test'))
        with self.assertRaises(AttributeError):
            self.structure.test


class TestDictRepr(unittest.TestCase):

    def setUp(self):
        """Create empty cell to work with."""
        self.structure = Structure()

    def test_attribute_in_dict(self):
        """Add attributes, remove attributes and test that they work"""
        self.structure.test = 101
        self.structure.nested = Structure()
        self.assertTrue('test' in self.structure.as_dict)
        self.assertTrue('nested' in self.structure.as_dict)
        self.assertIsInstance(self.structure.as_dict['nested'], dict)


if __name__ == '__main__':
    unittest.main(verbosity=3)
