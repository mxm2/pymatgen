#!/usr/bin/env python

"""
This module provides basic manipulation with siesta input and
output files.
"""


from __future__ import division

__author__ = "Artem Prihodko"
__version__ = "0.2"
__maintainer__ = "Artem Prihodko"
__email__ = "artemcpp@gmail.com"
__date__ = "Jun, 2014"

from pymatgen.core import Element, Molecule, Composition, Structure, Lattice
import re
import numpy as np


class SiestaInput(object):
	"""
    An object representing a Siesta input file.

    Args:
        structure: input structure.

        name: descriptive name of the system.

        spin_list: a list of spins for each element of the structure in the correct order.
    """

    # Pre-defined regex patterns
	params_pattern = re.compile('^\s*([\w.]+)\s+(.+)', re.M | re.I)
	block_pattern = re.compile('^\s*%block\s+(?P<n>[\w.]*)(.*)%endblock *(?P=n)',
									re.DOTALL | re.M | re.I)

	ScaleConstant = {'Ang': 1, 'Bohr': 0.5291772083}


	def __init__(self, structure, name = None, spin_list = None):
		self._structure = structure
		self._name = name
		self._spin_list = spin_list


	@property
	def structure(self):
		"""
		Returns structure associated with this SiestaInput
		"""
		return self._structure
	
	@property
	def name(self):
		"""
		Returns name of the system.
		"""
		return self._name

	@property
	def spin_list(self):
		"""
		Returns a list with spin values for each element.
		"""
		return self._spin_list


	@staticmethod
	def MakeTable(arr, ncolumns):
		"""
		Reshaping raw data in a table with n columns.
		"""
		if arr.size % ncolumns == 0:
			return arr.reshape(arr.size / ncolumns, ncolumns)
		else:
			return None


	@staticmethod
	def FromString(contents):
		"""
		Creates SiestaInput from a string.
		Args:
            contents: String representing an Siesta input file.

        Returns:
            SiestaInput object
		"""

		#TODO: Reading fdf ignoring - or _ like: LatticeConstant lattice_constant

		params = SiestaInput.params_pattern.findall(contents)
		blocks = SiestaInput.block_pattern.findall(contents)

		# Making Dictionary from raw list, deliting all \n entries
		block_dict = {blocks[i][0]: blocks[i][1].strip().replace('\n', ' ') for i in range(0, len(blocks))}
		param_dict = {params[i][0]: params[i][1].strip() for i in range(0, len(params))}

		if 'LatticeVectors' in block_dict.keys():
			lattice_list = np.array(block_dict['LatticeVectors'].strip().split()).astype(float)
			lattice = Lattice(lattice_list)
		elif 'LatticeParameters' in block_dict.keys():
			lattice_list = np.array(block_dict['LatticeParameters'].strip().split()).astype(float)
			lattice = Lattice.from_lengths_and_angles(lattice_list[:3], lattice_list[3:])


		# Getting Coordinates matrix
		if 'AtomicCoordinatesAndAtomicSpecies' in block_dict.keys():
			spec_data = block_dict['AtomicCoordinatesAndAtomicSpecies']		
			npspec_data = np.array(spec_data.split()).astype(float)
			npspec_table = SiestaInput.MakeTable(npspec_data, 4)
			npspec_table = npspec_table[np.lexsort(npspec_table.T)]	


		# Getting labels of Chemical elements
		chem_spec = np.array(block_dict['ChemicalSpeciesLabel'].split())
		chem_spec = SiestaInput.MakeTable(chem_spec, 3)

		elems = []
		for x in npspec_table:
			elems.append(chem_spec[x[3] - 1][2])


		# Checking Coordinates format and converting bohr to ang if necessery
		lattice_constant = param_dict['LatticeConstant'].strip().split()
		coord_format = param_dict['AtomicCoordinatesFormat']
		in_cartasian = False

		if coord_format:
			if coord_format in ('Ang', 'NotScaledCartesianAng'):
				in_cartasian = True;
			elif coord_format in ('Bohr', 'NotScaledCartesianBohr'):
				npspec_table[:, :-1] *= SiestaInput.ScaleConstant['Bohr']
			elif coord_format == 'ScaledCartesian':
				in_cartasian = True;
				npspec_table[:, :-1] *= SiestaInput.ScaleConstant[lattice_constant[1]] * float(lattice_constant[0])


		struct = Structure(lattice, elems, npspec_table[:, :-1], in_cartasian)

		#Parse initspin block:
		spin_list = None
		if 'DM.initspin' in block_dict.keys():
			spin_data = block_dict['DM.initspin']
			spin_list = spin_data.split()
			del spin_list[::2]
			spin_list = [float(spin) for spin in spin_list]


		return SiestaInput(struct, param_dict["SystemName"], spin_list)

	@staticmethod
	def from_file(filename):
		"""
        Creates SiestaInput from a file.

        Args:
            filename: Siesta input filename

        Returns:
            SiestaInput object
        """
		with open(filename) as f:
			return SiestaInput.FromString(f.read())


class SiestaOutput(object):
	"""
	An object representing Siesta output file.

	Args:
		filename: path to Siesta output file.
	"""

	def __init__(self, filename):
		self._filename = filename
		self.parse(filename)


	def parse(self, filename):
		"""
		Parsing function for standard Siesta output.

		Args:
			filename: path to Siesta output file.
		"""
		def _find(pattern, str):
			match = pattern.search(str)
			if match:
				return match.group(1)
			else:
				return None

		ver_pattern = re.compile("Siesta Version\s*:\s*(.*)", re.I)
		arc_pattern = re.compile("Architecture\s*:\s*(.*)", re.I | re.M)
		nspecies_pattern = re.compile("Number of Atomic Species\s*=\s*(\d+)", re.I | re.M)
		stype_pattern = re.compile("System type\s*=\s*(\w+)", re.I | re.M)
		num_kpts_pattern = re.compile("Number of k-points\s*=\s*(\d+)", re.I | re.M)
		mesh_cutoff_pattern = re.compile("Mesh cutoff \(required, used\)\s*=\s*(.+)Ry", re.I | re.M)

		with open(filename) as f:
			contents = f.read()

			self._ver = _find(ver_pattern, contents)
			self._arc = _find(arc_pattern, contents)
			self._nspecies = int(_find(nspecies_pattern, contents))
			self._sys_type = _find(stype_pattern, contents)
			self._num_kpts = int(_find(num_kpts_pattern, contents))
			self._mesh_cutoff = [float(mk) for mk in _find(mesh_cutoff_pattern, contents).split()]



	@property
	def ver(self):
	    return self._ver

	@property
	def arc(self):
	    return self._arc

	@property
	def nspecies(self):
	    return self._nspecies

	@property
	def sys_type(self):
	    return self._sys_type

	@property
	def num_kpts(self):
	    return self._num_kpts

	@property
	def mesh_cutoff(self):
	    return self._mesh_cutoff