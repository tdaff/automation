"""
Unit testing based on set of test cases.

"""

from .. import faps
import unittest

class TestStructureReaders(unittest.TestCase):

    def setUp(self):
        self.structure = faps.Structure('test')

    def test_pdb(self):
        # make sure the shuffled sequence does not lose any elements
        self.structure.from_pdb('test.pdb')
        random.shuffle(self.seq)
        self.seq.sort()
        self.assertEqual(self.seq, range(10))

        # should raise an exception for an immutable sequence
        self.assertRaises(TypeError, random.shuffle, (1,2,3))

    def test_choice(self):
        element = random.choice(self.seq)
        self.assertTrue(element in self.seq)

    def test_sample(self):
        with self.assertRaises(ValueError):
            random.sample(self.seq, 20)
        for element in random.sample(self.seq, 5):
            self.assertTrue(element in self.seq)

if __name__ == '__main__':
    unittest.main()
