#!/usr/bin/python

import os
import sys
import urllib2
import pymatgen.matproj.rest as rester
import pymatgen.io.smartio as sio
import pymatgen.io.vaspio_set

#Allows to convert .cif files into siesta input file using materialsproject vasp settings.
#Early test version. Going to be modified and commented properly.
#Usage: Place ciffiles to any folder. 'cd' this folder and run gen_siesta_input.py *.cif

API_KEY = "UvSgztFQ7CcO6fb5"
PSF_LINK = "http://departments.icmab.es/leem/siesta/Databases/Pseudopotentials/Pseudos_GGA_Abinit/"


def get_full_vasp_set():
	mp_api = rester.MPRester(API_KEY)
	return mp_api.get_vasp_input_set()

def get_structure_from_file(filename):
	return sio.read_structure(filename)

def get_psf_files(elem_list, fdir):
	for elem in elem_list:
		try:
			psf = urllib2.urlopen(PSF_LINK + elem + '_html/' + elem + '.psf')
		except:
			print "Warning: Couldn't find psf for '%s'" % elem
			continue
		with open(fdir + '/' + elem + '.psf', 'w') as f:
			f.write(psf.read())

def MakeBlock(block_name, content):
	block = "%block " + block_name + "\n"
	for el in content:
		block += " ".join(list(el)) + "\n"
	block += "%endblock " + block_name + "\n"
	return block

def write_to_vasp(structure, vasp, directory):
	vasp.get_incar(structure).write_file(directory + '/INCAR')
	vasp.get_poscar(structure).write_file(directory + '/POSCAR')
	vasp.get_kpoints(structure).write_file(directory + '/KPOINTS')

def write_to_fdf(structure, vasp, fdf_filename):
	incar = vasp.get_incar(structure)
	kpoints = vasp.get_kpoints(structure).to_dict
	kgrid = kpoints['kpoints'][0]
	shift = kpoints['usershift']
	struct_dict = structure.to_dict

	elems = dict(zip(structure.symbol_set, range(1, len(structure.symbol_set) + 1)))
	anums = structure.atomic_numbers
	anums = [v for i,v in enumerate(anums) if v not in anums[:i]]
	ChemSpec = [(str(n), str(anum), elem) for elem, n, anum in zip(elems.iterkeys(), elems.itervalues(), anums)]

	AtomicCoord = [["%.5f %.5f %.5f    " % tuple(site['abc']) + str(elems[site["label"]])] for site in struct_dict['sites']]


	with open(fdf_filename, 'w') as fdf_file:
		fdf_file.write("SystemName %s\n" % structure.formula)
		fdf_file.write("SystemLabel %s\n\n" % structure.formula)
		fdf_file.write("NumberOfAtoms %s\n" % structure.num_sites)
		fdf_file.write("NumberOfSpecies %s\n\n" % structure.ntypesp)
		fdf_file.write(MakeBlock("ChemicalSpeciesLabel", ChemSpec))
		fdf_file.write("LatticeConstant 1 Ang\n\n")
		fdf_file.write(MakeBlock("LatticeParameters", [["%(a).4f %(b).4f %(c).4f    %(alpha).4f %(beta).4f %(gamma).4f" % struct_dict['lattice']]]))
		fdf_file.write("AtomicCoordinatesFormat  Fractional\n")
		fdf_file.write(MakeBlock("AtomicCoordinatesAndAtomicSpecies", AtomicCoord))
		fdf_file.write("ElectronicTemperature %(SIGMA)s eV\n\n" % incar)
		
		if int(incar["ISMEAR"]) in range(-5, 0):
			fdf_file.write("OccupationFunction FD\n\n")
		elif int(incar["ISMEAR"]) >= 1:
			fdf_file.write("OccupationFunction MP\n")
			fdf_file.write("OccupationMPOrder %(ISMEAR)s\n\n" % incar)	

		if set(incar['MAGMOM']) != {0.6}:
			fdf_file.write("SpinPolarised .true.\n")
			fdf_file.write(MakeBlock("DM.initspin", [["%i    %.1f" % tuple([i + 1, l])] for i, l in enumerate(incar['MAGMOM'])]))
		fdf_file.write('\n')

		# if kpoints['generation_style'] == 'Gamma':
			# shift = [i - 0.5 for i in shift]
		fdf_file.write(MakeBlock("kgrid_Monkhorst_pack", [["%i 0 0 " % kgrid[0] + str(shift[0]) + "\n" +\
		 "0 %i 0 " % kgrid[1] + str(shift[1]) + "\n" +\
		 "0 0 %i " % kgrid[2] + str(shift[2])]]))
		fdf_file.write('\n')

		fdf_file.write("MD.NumCGsteps       100\n")
		fdf_file.write("PAO.EnergyShift     65 meV\n")
		fdf_file.write("MeshCutoff          300 Ry\n")
		fdf_file.write("MaxSCFIterations    1000\n\n")
		fdf_file.write("LongOutput .true.\n")
		fdf_file.write("WriteForces .true.")



if len(sys.argv) < 2:
	print "Error: First argument must be the name of cif file."
	sys.exit()

for cif_filename in sys.argv[1:]:
	if not os.path.exists(cif_filename):
		print "Error: File '%s' does not exist." % sys.argv[1]
		sys.exit()

	vasp = get_full_vasp_set()
	structure = get_structure_from_file(cif_filename)

	dirname = cif_filename.split('/')[-1]
	dirname = dirname.split('.')[0]
	if not os.path.exists(dirname):
		os.makedirs(dirname)

	# write_to_vasp(structure, vasp, dirname)
	write_to_fdf(structure, vasp, dirname + '/t4x.fdf')
	get_psf_files(structure.symbol_set, dirname)

	print "Done: " + dirname
