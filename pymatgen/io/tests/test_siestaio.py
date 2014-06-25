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
from pymatgen.io.siestaio import SiestaInput, SiestaOutput

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

class SiestaInputTest(unittest.TestCase):
	#TODO: Set correct Fe_output attributes

	def SetUp(self):
		pass

	def test_from_file(self):
		Fe = SiestaInput.from_file(os.path.join(test_dir, 'Fe.fdf'))
		Na2O = SiestaInput.from_file(os.path.join(test_dir, 'Na2O1.fdf'))

		self.assertAlmostEqual(Fe.structure.volume, 256.307, places = 3)
		self.assertEqual(Fe.name, 'Fe1')
		self.assertEqual(Fe.spin_list, [5.0])

		self.assertAlmostEqual(Na2O.structure.volume, 43.789, places = 3)
		self.assertEqual(Na2O.name, 'Na2O1')
		self.assertEqual(Na2O.spin_list, [0.6, 0.6, 0.6])

class SiestaOutputTest(unittest.TestCase):
	def SetUp(self):
		pass

	def test_props(self):
		feout = SiestaOutput(os.path.join(test_dir, 'Fe_output'))

		self.assertEqual(feout.ver, 'siesta-3.2-pl-4')
		self.assertEqual(feout.arc, 'i686-pc-linux-gnu--unknown')
		self.assertEqual(feout.nspecies, 1)
		self.assertEqual(feout.sys_type, 'bulk')
		self.assertEqual(feout.num_kpts, 560)
		self.assertEqual(feout.mesh_cutoff[0], 300)
		self.assertEqual(feout.mesh_cutoff[1], 400.750) #in Ry





if __name__ == "__main__":
	unittest.main()