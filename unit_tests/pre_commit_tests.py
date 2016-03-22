import cPickle as pickle
from merlin import *

def test_xyz_cml():
	molecule = utils.Molecule('acetone')
	atoms = files.read_xyz('PbCl_24')
	for filetype in ['xyz', 'cml']:
		if os.path.isfile('out.'+filetype): os.remove('out.'+filetype)
		if filetype=='xyz': files.write_xyz(atoms)
		else: files.write_cml(molecule)
		if open('out.'+filetype).read() != open('target_out.'+filetype).read():
			raise Exception(filetype+' filetype test failed')

def test_orca():
	result = orca.read('PbCl2_0_vac')
	target_result = pickle.load(open('orca.pickle'))
	if result.energy != target_result.energy:
		raise Exception('Wrong orca energy: %f vs %f' % (result.energy, target_result.energy))

def test_g09():
	result = g09.read('PbCl2_0_vac')
	target_result = pickle.load(open('g09.pickle'))
	if result.energy != target_result.energy:
		raise Exception('Wrong g09 energy: %f vs %f' % (result.energy, target_result.energy))

def test_lammps():
	molecule = utils.Molecule('acetone')
	system = utils.System(box_size=[30.0,30.0,30.0], name='test')
	system.add(molecule)
	os.chdir('lammps')
	files.write_lammps_data(system)
	os.chdir('..')

def test_files():
	os.chdir('unit_tests/test_files')
	test_xyz_cml()
	test_orca()
	test_g09()
	test_lammps()
	os.chdir('..')

test_files()

print "All pre-commit tests succeeded"


