import cPickle as pickle
import re, shutil
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
	orca.job('test', 'HF-3c', previous='H2').wait()
	result = orca.read('test')
	target_result = orca.read('H2')
	if result.energy != target_result.energy:
		raise Exception('Wrong orca energy: %f vs %f' % (result.energy, target_result.energy))

def test_g09():
	result = g09.read('PbCl2_0_vac')
	target_result = pickle.load(open('g09.pickle'))
	if result.energy != target_result.energy:
		raise Exception('Wrong g09 energy: %f vs %f' % (result.energy, target_result.energy))

def test_lammps():
	molecule = utils.Molecule('acetone')
	system = utils.System(box_size=[20.0,20.0,20.0], name='test')
	system.add(molecule)
	files.packmol(system, (molecule,), molecule_ratio=(1,), density=0.1)
	if os.path.isdir("lammps"):
		shutil.rmtree('lammps')
	os.mkdir('lammps')
	os.chdir('lammps')
	files.write_lammps_data(system, pair_coeffs_included=True)
	output = open(system.name+'.in', 'w')
	output.write('''	units real
atom_style full
pair_style lj/cut/coul/dsf 0.05 2.5 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5
boundary p p p
read_data	'''+system.name+'''.data
thermo_style multi
thermo 0
velocity all create 300.0 1 rot yes dist gaussian
fix motion all nvt temp 300.0 300.0 100.0
run 10''')
	output.close()
	os.system('/fs/home/hch54/lammps/lammps-7Dec15/src/lmp_serial -in %s.in -log %s.log > /dev/null' % (system.name,system.name))
	logfile = open(system.name+'.log').read()
	energies = re.findall('TotEng += +(\S+)', logfile)
	target_energies = ['8.4933', '10.3644']
	if energies != target_energies:
		raise Exception('LAMMPS test failed: %s != %s' % (str(energies), str(target_energies)) )
	os.chdir('..')

def test_utils():
	molecule = utils.Molecule('acetone')
	molecule2 = utils.Molecule('acetone')
	assert molecule.equals(molecule2), "Molecule.equals() has failed tests."
	system = utils.System(box_size=[20.0,20.0,20.0], name='test')
	system.add(molecule)
	system2 = utils.System(box_size=[20.0,20.0,20.0], name='test')
	system2.add(molecule2)
	assert system.equals(system2), "System.equals() has failed tests."
	assert system.Contains(molecule), "System.Contains() has failed tests."
	system.Remove(molecule)
	molecule.translate([1.0,0.0,0.0])

	assert not molecule.equals(molecule2), "Molecule.equals() has failed tests."
	assert not system.Contains(molecule), "System.Contains() has failed tests."
	assert not system2.Contains(molecule), "System.Contains() has failed tests."
	system.add(molecule)
	system.Remove(molecule)
	
	system3 = utils.System(box_size=[20.0,20.0,20.0], name = 'test')
	assert system.equals(system3), "System.Remove() has failed tests."
	

def test_files():
	os.chdir('unit_tests/test_files')
	test_xyz_cml()
	test_orca()
	test_g09()
	test_lammps()
	test_utils()
	os.chdir('..')

test_files()

print "All pre-commit tests succeeded"


