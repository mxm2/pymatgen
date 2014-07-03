from pymatgen.core import Element, Molecule, Composition, Structure, Lattice
import re
import numpy as np

class SiestaInput(object):
	# def __init(self, mol, charge = None, title = None, 
	# 		   input_parameters = None):
	# 	self._mol = mol
	# 	self.charge = charge if charge in not None else mol.charge

	def FromString(self, contents):
		ScaleConstant = {'Ang': 1, 'Bohr': 0.5291772083}

		params_pattern = re.compile('^([A-Za-z]*) (.*)\n', re.M)
		block_pattern = re.compile('^%block (?P<n>[A-Za-z_]*)(.*)%endblock *(?P=n)',
									re.DOTALL | re.M)
		params = params_pattern.findall(contents)
		blocks = block_pattern.findall(contents)


		# Making Dictionary from raw list, deliting all \n entries
		block_dict = {blocks[i][0]: blocks[i][1].strip().replace('\n', ' ') for i in range(0, len(blocks))}
		param_dict = {params[i][0]: params[i][1].strip() for i in range(0, len(params))}

		# print "Full FDF Dictionary:\n", block_dict

		if 'LatticeVectors' in block_dict.keys():
			lattice_list = np.array(block_dict['LatticeVectors'].strip().split()).astype(float)
			lattice = Lattice(lattice_list)
		elif 'LatticeParameters' in block_dict.keys():
			lattice_list = np.array(block_dict['LatticeParameters'].strip().split()).astype(float)
			lattice = Lattice.from_lengths_and_angles(lattice_list[:3], lattice_list[3:])
		print lattice_list[:3], '\n', lattice_list[3:]

		# print "Lattice Matrix:\n" , lattice.matrix
		# print "Siesta input parameters:\n", param_dict

		# Getting Coordinates matrix
		if 'AtomicCoordinatesAndAtomicSpecies' in block_dict.keys():
			Coords = block_dict['AtomicCoordinatesAndAtomicSpecies']		
			npCoord = np.array(Coords.split()).astype(float)
			if npCoord.size % 4 == 0:
				npCoord = npCoord.reshape(npCoord.size / 4, 4)
				# Sorting table by element number
				npCoord = npCoord[np.lexsort(npCoord.T)]

		# Getting labels of Chemical elements
		ChemSpec = np.array(block_dict['ChemicalSpeciesLabel'].split())
		if ChemSpec.size % 3 == 0:
			ChemSpec = ChemSpec.reshape(ChemSpec.size / 3, 3)

		Elems = []
		for x in npCoord:
			Elems.append(ChemSpec[x[3] - 1][2])

		# Checking Coordinates format and converting bohr to ang if necessery
		LatticeConstant = param_dict['LatticeConstant'].strip().split()
		inCartesian = False
		CoordFormat = param_dict['AtomicCoordinatesFormat']
		if CoordFormat:
			if CoordFormat in ('Ang', 'NotScaledCartesianAng'):
				inCartesian = True;
			elif CoordFormat in ('Bohr', 'NotScaledCartesianBohr'):
				inCartesian = True;
				npCoord[:, :-1] *= ScaleConstant['Bohr']
			elif CoordFormat == 'ScaledCartesian':
				inCartesian = True;
				npCoord[:, :-1] *= ScaleConstant[LatticeConstant[1]] * float(LatticeConstant[0])

		struct = Structure(lattice, Elems, npCoord[:, :-1], inCartesian)
		print '---------------\n', struct.to_dict, '\n---------------\n'

		return contents


	def FromFile(self, filename):
		with open(filename) as f:
			return self.FromString(f.read())


siesta = SiestaInput()

for f in ('/home/artem/nano/Fe_isolated_spin_5/Fe_isolated/t4x.fdf', '/home/artem/nano/Na.fdf', '/home/artem/nano/Na2', '/home/artem/nano/Na3', '/home/artem/nano/Na4'):
	siesta.FromFile(f)
