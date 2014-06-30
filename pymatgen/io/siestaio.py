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

    # Pre-defined regex patts
	params_patt = re.compile('^\s*([\w.]+)\s+(.+)', re.M | re.I)
	block_patt = re.compile('^\s*%block\s+(?P<n>[\w.]*)(.*)%endblock *(?P=n)',
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

		params = SiestaInput.params_patt.findall(contents)
		blocks = SiestaInput.block_patt.findall(contents)

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

	Attributes:

		.. attribute:: ver

			Version of Siesta used for run.

		.. attribute:: arc

			Architecture. 

		.. attribute:: nspecies

			Number of Species.

		.. attribute:: sys_type

			Type of th system (bulk, surface, single atom)

		.. attribute:: num_atoms

			Number of atomsin the run.

		.. attribute:: num_kpts

			Number of k-points in run

		.. attribute:: mesh_cutoff

			List of 2 elements:
		 		First: required value of MeshCutoff (Ry)
		 		Second: used value (Ry)

		.. attribute:: total_energy

			Total energy of the system in eV

		.. attribute:: charge

			Charge for structure.

		.. attribute:: num_basis_func

			Number of basis functions in the run.

		.. attribute:: pressure_solid

			Static pressure for solid in Ry/Bohr**3.

		.. attribute:: pressure_mol

			Static pressure for molecule Ry/Bohr**3.

		.. attribute:: cell_vol

			Cell volume in Ang**3.

		.. attribute:: atomic_force

			Max atomic force in eV/Ang.
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
		def find_by_pattern(pattern, str):
			match = pattern.search(str)
			if match:
				return match.group(1)
			else:
				return None

		def cut(content, begin_str, end_str):
			pos_begin = content.find(begin_str)
			pos_end = content.find(end_str, pos_begin)
			if pos_begin != -1 and pos_end > pos_begin:
				return content[pos_begin:pos_end]
			else:
				return ''

		ver_patt = re.compile("Siesta Version\s*:\s*(.*)", re.I)
		arc_patt = re.compile("Architecture\s*:\s*(.*)", re.I | re.M)
		nspecies_patt = re.compile("Number of Atomic Species\s*=\s*(\d+)", re.I | re.M)
		stype_patt = re.compile("System type\s*=\s*(\w+)", re.I | re.M)
		natom_patt = re.compile("superc: Number of atoms, orbitals, and projectors\s*:\s*(\d+)\s+", re.I | re.M)
		num_kpts_patt = re.compile("Number of k-points\s*=\s*(\d+)", re.I | re.M)
		mesh_cutoff_patt = re.compile("Mesh cutoff \(required, used\)\s*=\s*(.+)Ry", re.I | re.M)
		f_energy_patt = re.compile("Final energy (eV):", re.I | re.M)
		qtot_patt = re.compile("Qtot\s*=\s*([\d.-]+)")
		basisfun_patt = re.compile("nzeta=(\d+)\s+polorb=(\d+)")
		press_patt = re.compile("([\d.]+)\s*([\d.]+)\s*Ry/Bohr\*\*3")
		cell_vol_patt = re.compile("Cell volume\s*=\s*([\d.]+)")
		maxforce_patt = re.compile("Max\s*([\d.-]+)")

		s_line_patt = re.compile(r"siesta:\s*(.+)\s*=\s*([\d.-]+)\s*")

		with open(filename) as f:
			content = f.read()

			self.ver = find_by_pattern(ver_patt, content)
			self.arc = find_by_pattern(arc_patt, content)
			self.nspecies = int(find_by_pattern(nspecies_patt, content))
			self.sys_type = find_by_pattern(stype_patt, content)
			self.num_atoms = int(find_by_pattern(natom_patt, content))
			self.num_kpts = int(find_by_pattern(num_kpts_patt, content))
			self.mesh_cutoff = [float(mk) for mk in find_by_pattern(mesh_cutoff_patt, content).split()]

			self.total_energy = False
			self.charge = False
			self.pressure_solid = False
			self.pressure_mol = False
			self.num_basis_func = -1
			self.cell_vol = -1

			#Parse Final Energy
			fe_block = cut(content, 'Final energy (eV):', 'siesta: Stress tensor')
			en_list = s_line_patt.findall(fe_block)
			if en_list:
				final_energy = {}
				for en_type, value in en_list:
					final_energy[en_type.strip()] = float(value)
				if 'Total' in final_energy:
					self.total_energy = final_energy['Total']
			
			#charge (in mulliken: spin up + spin down)
			up_down_charge = qtot_patt.findall(content)
			if len(up_down_charge) == 2:
				self.charge = float(up_down_charge[0]) + float(up_down_charge[1])

			#number of basis functions
			basis_block = cut(content, '<basis_specs>', '</basis_specs>')
			basis_data = basisfun_patt.findall(basis_block)
			if basis_data:
				nfunc = 0
				for nzeta, polarb in basis_data:
					nfunc = nfunc + int(nzeta) + int(polarb)
			if nfunc:
				self.num_basis_func = nfunc * self.num_atoms

			#pressure in Ry/Bohr**3
			pressure = {}
			press_block = cut(content, 'Pressure (static):', '(Free)E+ p_basis*V_orbitals')
			press_data = press_patt.search(press_block)
			if press_data:
				self.pressure_solid = float(press_data.group(1))
				self.pressure_mol = float(press_data.group(2))

			#cell volume in Ang**3
			cell_vol = cell_vol_patt.search(content)
			if cell_vol:
				self.cell_vol = float(cell_vol.group(1))

			#max atomic force in eV/Ang
			aforces_block = cut(content, 'Atomic forces (eV/Ang):', 'Stress-tensor-Voigt')
			aforces_data = maxforce_patt.findall(aforces_block)
			if aforces_data:
				self.atomic_force = float(aforces_data[-1])
