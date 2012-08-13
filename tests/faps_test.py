#!usr/bin/env python

"""
Unit testing based on set of test cases.

"""

import sys
import unittest

# This is where faps gets imported
sys.path.append('..')
import config
import faps

class TestConfigReader(unittest.TestCase):

    def setUp(self):
        self.options = config.Options(job_name='test')

    def test_options(self):
        self.assertEqual(self.options.get('job_name'), 'test')


class TestStructureReaders(unittest.TestCase):

    def setUp(self):
        self.structure = faps.Structure('test')

    def test_pdb(self):
        # make sure the shuffled sequence does not lose any elements
        #self.structure.from_pdb('test.pdb')
        pass


class TestCellCrystalSystems(unittest.TestCase):

    def setUp(self):
        self.cell = faps.Cell()

    def test_crystal_perception(self):
        self.cell.params = (10, 10, 10, 90, 90, 90)
        self.assertEqual(self.cell.crystal_system, 'cubic')
        self.cell.params = (10, 10, 12, 90, 90, 90)
        self.assertEqual(self.cell.crystal_system, 'tetragonal')
        self.cell.params = (10, 11, 10, 90, 90, 90)
        self.assertEqual(self.cell.crystal_system, 'tetragonal')
        self.cell.params = (13, 10, 10, 90, 90, 90)
        self.assertEqual(self.cell.crystal_system, 'tetragonal')
        self.cell.params = (10, 11, 12, 90, 90, 90)
        self.assertEqual(self.cell.crystal_system, 'orthorhombic')
        self.cell.params = (10, 11, 12, 90, 70, 90)
        self.assertEqual(self.cell.crystal_system, 'monoclinic')
        self.cell.params = (10, 11, 12, 90, 90, 70)
        self.assertEqual(self.cell.crystal_system, 'monoclinic')
        self.cell.params = (10, 11, 12, 70, 90, 90)
        self.assertEqual(self.cell.crystal_system, 'monoclinic')
        self.cell.params = (10, 11, 12, 60, 70, 80)
        self.assertEqual(self.cell.crystal_system, 'triclinic')
        self.cell.params = (10, 10, 12, 90, 90, 120)
        self.assertEqual(self.cell.crystal_system, 'hexagonal')
        self.cell.params = (10, 11, 10, 90, 120, 90)
        self.assertEqual(self.cell.crystal_system, 'hexagonal')
        self.cell.params = (13, 10, 10, 120, 90, 90)
        self.assertEqual(self.cell.crystal_system, 'hexagonal')
        self.cell.params = (10, 10, 10, 67, 67, 67)
        self.assertEqual(self.cell.crystal_system, 'trigonal')

if __name__ == '__main__':
    # verbosity argument only added in 2.7
    if sys.version_info[:2] < (2, 7):
        unittest.main()
    else:
        unittest.main(verbosity=2)
