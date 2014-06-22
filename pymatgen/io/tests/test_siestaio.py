#!/usr/bin/env python

"""
Created on Jun 15, 2014
"""


from __future__ import division

__author__ = "Artem Prihodko"
__version__ = "0.1"
__maintainer__ = "Artem Prihodko"
__email__ = "artemcpp@gmail.com"
__date__ = "Jun, 2014"


import unittest
import os

from pymatgen.core import Element, Molecule, Composition, Structure, Lattice
from pymatgen.io.siestaio import SiestaInput

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "fdf")


class SiestaInputTest(unittest.TestCase):
	def SetUp(self):
		pass

	def test_from_file(self):
		Fe = SiestaInput.from_file(os.path.join(test_dir, 'Fe.fdf'))
		Na2O = SiestaInput.from_file(os.path.join(test_dir, 'Na2O1.fdf'))

		self.assertAlmostEqual(Fe.structure.volume, 256.307, places = 3)
		self.assertEqual(Fe.label, 'Fe1')
		self.assertEqual(Fe.spin_list, [5.0])

		self.assertAlmostEqual(Na2O.structure.volume, 43.789, places = 3)
		self.assertEqual(Na2O.label, 'Na2O1')
		self.assertEqual(Na2O.spin_list, [0.6, 0.6, 0.6])



if __name__ == "__main__":
	unittest.main()