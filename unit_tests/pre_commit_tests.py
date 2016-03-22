from merlin import *

def test_files():
	os.chdir('unit_tests/test_files')
	molecule = utils.Molecule('acetone')
	atoms = files.read_xyz('PbCl_24')
	for filetype in ['xyz', 'cml']:
		if os.path.isfile('out.'+filetype): os.remove('out.'+filetype)
		if filetype=='xyz': files.write_xyz(atoms)
		else: files.write_cml(molecule)
		if open('out.'+filetype).read() != open('target_out.'+filetype).read():
			raise Exception(filetype+' filetype test failed')
	os.chdir('..')

test_files()

print "All pre-commit tests succeeded"


