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
		fe_out = SiestaOutput(os.path.join(test_dir, 'Fe_output'))

		self.assertEqual(fe_out.ver, 'siesta-3.2-pl-4')
		self.assertEqual(fe_out.arc, 'i686-pc-linux-gnu--unknown')
		self.assertEqual(fe_out.nspecies, 1)
		self.assertEqual(fe_out.sys_type, 'bulk')
		self.assertEqual(fe_out.num_atoms, 1000)
		self.assertEqual(fe_out.num_kpts, 560)
		self.assertEqual(fe_out.mesh_cutoff[0], 300)
		self.assertEqual(fe_out.mesh_cutoff[1], 400.750) #in Ry
		self.assertEqual(fe_out.total_energy, -567.924870)
		self.assertEqual(fe_out.charge, 8)
		# # self.assertEqual(fe_out.spin_multiplicity, 8)
		self.assertEqual(fe_out.num_basis_func, 5000)
		self.assertEqual(fe_out.pressure_solid, 0.00209287)
		self.assertEqual(fe_out.pressure_mol, 0.00209287)
		self.assertEqual(fe_out.cell_vol, 11.196570)
		self.assertEqual(fe_out.atomic_force, 0.000005)


		graphene_out = SiestaOutput(os.path.join(test_dir, 'graphene_Sc_output'))

		self.assertEqual(graphene_out.ver, 'siesta-3.2-pl-4')
		self.assertEqual(graphene_out.arc, 'x86_64-unknown-linux-gnu--unknown')
		self.assertEqual(graphene_out.nspecies, 2)
		self.assertEqual(graphene_out.sys_type, 'slab')
		# self.assertEqual(graphene_out.num_atoms, 864)
		self.assertEqual(graphene_out.num_kpts, 1152)
		self.assertEqual(graphene_out.mesh_cutoff[0], 300.000)
		self.assertEqual(graphene_out.mesh_cutoff[1], 314.455) #in Ry
		self.assertEqual(graphene_out.total_energy, -1138.731881)
		# self.assertEqual(graphene_out.charge, 31) #! return False
		# # self.assertEqgraphenel(fe_out.spin_multiplicity, 8)
		# self.assertEqual(graphene_out.num_basis_func, 5*864)
		self.assertEqual(graphene_out.pressure_solid, 0.00000210)
		self.assertEqual(graphene_out.pressure_mol, 0.00000200)
		self.assertEqual(graphene_out.cell_vol, 377.651841)
		self.assertEqual(graphene_out.atomic_force, 0.009650)


if __name__ == "__main__":
	unittest.main()