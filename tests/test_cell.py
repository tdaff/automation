#!usr/bin/env python

"""
Testing Cell functionality

"""

import sys
import unittest

# This is where faps gets imported
sys.path.append('..')
from faps.core.cell import Cell

class TestInitialiseCell(unittest.TestCase):

    def test_initialise_cell(self):
        """Just make a new cell and test the defaults."""
        new_cell = Cell()
        self.assertEqual(new_cell.cell.tolist(),
                         [[1.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0],
                          [0.0, 0.0, 1.0]])
        self.assertEqual(new_cell.params, (1.0, 1.0, 1.0, 90, 90, 90))


class TestCellCrystalSystems(unittest.TestCase):

    def setUp(self):
        # cell to manipulate
        self.cell = Cell()

    def test_crystal_perception(self):
        # Make sure that DL_POLY imcons work
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
    unittest.main(verbosity=3)
