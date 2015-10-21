from merlin import *
from subprocess import Popen

# A function to run an Orca DFT Simulation
def job(run_name, route, atoms=[], extra_section='', grad=False, queue=None, procs=1, charge_and_multiplicity='0 1', title='', blurb=None, force=False, previous=None, neb=[False,None,None,None], lint=False):
	# Generate the orca input file
	os.system('mkdir -p orca/%s' % run_name)
	os.chdir('orca/%s' % run_name)

	# If running on system with more than one core, tell orca
	if procs > 1:
		route += ' PAL%d' % procs

	if grad:
		extra_section = extra_section.strip() + '''\n%method
 RunTyp Gradient
 end'''

	# Get input for orca formatted correctly
	inp = route.strip()+'\n'+extra_section.strip()+'\n*xyz '+charge_and_multiplicity+'\n'
	for a in atoms:
		inp += '%s %f %f %f\n' % (a.element, a.x, a.y, a.z)
	inp += '*\n'

	# Write the orca file
	f = open(run_name+'.orca', 'w')
	f.write(inp)
	f.close()

	# Run the simulation
	if queue is None:
		process_handle = Popen('/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca %s.orca > %s.out' % (run_name, run_name), shell=True)
	else:
		NBS = '''#!/bin/bash
##NBS-name: "%s"
##NBS-nproc: %d
##NBS-queue: "%s"

export PATH=/fs/europa/g_pc/ompi_1_6_5/bin:/fs/europa/g_pc/orca_3_0_3_linux_x86-64:$PATH
export LD_LIBRARY_PATH=/fs/europa/g_pc/ompi_1_6_5/lib:$LD_LIBRARY_PATH

/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca %s.orca >> %s.out 2>&1
''' % (run_name, procs, queue, os.getcwd()+'/'+run_name, os.getcwd()+'/'+run_name)
		f = open(run_name+'.nbs', 'w')
		f.write(NBS)
		f.close()
		os.system('jsub %s.nbs' % run_name)

	# Return to the appropriate directory
	os.chdir('../..')

	if queue is None:
		return process_handle
	else:
		return utils.Job(run_name)

# A function to parse and Orca DFT output file (assumes by default a .out file format)
def parse_atoms(input_file, get_atoms=True, get_energy=True, get_charges=False, charge_type='MULLIKEN', get_time=False, check_convergence=False, parse_all=False):
	data = open('orca/%s/%s.out' % (input_file,input_file),'r').read()

	# Get all the positions
	if get_atoms:
		hold, atoms = data, []
		s = 'CARTESIAN COORDINATES (ANGSTROEM)'
		while hold.find(s) != -1:
			hold = hold[hold.find(s)+len(s):]
			tmp = hold[:hold.find('\n\n')].split('\n')[2:]
			tmp_atoms = []
			for a in tmp:
				a = a.split()
				tmp_atoms.append(utils.Atom(a[0],float(a[1]),float(a[2]),float(a[3])))
			atoms.append(tmp_atoms)
		if not parse_all: atoms = atoms[-1]

	# Get all the energies
	if get_energy:
		hold, energies = data, []
		s = 'FINAL SINGLE POINT ENERGY'
		while hold.find(s) != -1:
			hold = hold[hold.find(s):]
			tmp = hold[:hold.find('\n')].split()[-1]
			energies.append(float(tmp))
			hold = hold[hold.find('\n'):]
		if not parse_all: energies = energies[-1]

	# Get charges
	if get_charges:
		hold, charges = data, []
		if charge_type.strip().upper() not in ['MULLIKEN','LOEWDIN']:
			print("Error - Requested %s charge, but not a valid charge to request.  Choose from MULLIKEN or LOEWDIN." % charge_type.strip().upper())
			sys.exit()
		s = charge_type.strip().upper()+' ATOMIC CHARGES'

		hold = hold[hold.rfind(s):]
		b = hold[:hold.find('\n\n')].split('\n')[2:-1]
		for a in b:
			a = a.split()
			charges.append([a[1],float(a[3])])
	
	# Get Simulation Times (note, if not parse_all, you get total simulation time)
	if get_time:
		hold, times = data, []
		s = 'Total SCF time'
		while hold.find(s) != -1:
			hold = hold[hold.find(s):]
			tmp = hold[:hold.find('\n')].split()
			times.append(float(tmp[3])*3600*24 + float(tmp[5])*3600 + float(tmp[7])*60 + float(tmp[9]))
			hold = hold[hold.find('\n'):]
		if not parse_all: times = sum(times)

	if check_convergence:
		hold, convergence = data, []
		s = 'Geometry convergence'
		while hold.find(s) != -1:
			hold = hold[hold.find(s)+len(s):]
			tmp = hold[:hold.find('Max(Bonds)')].split('\n')[3:-2]
			tmp_convergence = []
			for a in tmp:
				a = a.split()
				tmp_convergence.append([' '.join(a[:2]),float(a[2]),float(a[3]), a[4]])
			convergence.append(tmp_convergence)
		if not parse_all: convergence = convergence[-1]

	results = []
	if get_atoms: results.append(atoms)
	if get_energy: results.append(energies)
	if get_charges: results.append(charges)
	if get_time: results.append(times)
	if check_convergence: results.append(convergence)
	return results

# A function to parse orca.engrad files
def engrad_read(input_file):
	if not os.path.isfile('orca/%s/%s.orca.engrad' % (input_file,input_file)):
		print("Currently in %s, searching for orca/%s/%s.orca.engrad" % (os.getcwd(),input_file,input_file))
		print("No engrad file exists. Please run simulation with grad=True.")
		sys.exit()
	data = open('orca/%s/%s.orca.engrad' % (input_file,input_file),'r').read().split('\n')
	count, grad, atoms = 0, [], []
	i = -1
	while i < len(data):
		i += 1
		line = data[i]
		if len(line) < 1: continue
		if line.strip()[0] == '#': continue
		if count == 0:
			num_atoms = int(line.strip())
			count += 1
		elif count == 1:
			# Get energy
			energy = float(line.strip())
			count += 1
		elif count == 2:
			# Get gradient
			for j in range(num_atoms):
				for k in range(3):
					grad.append(float(data[i+k].strip()))
				i += 3
			count += 1
		elif count == 3:
			# Get atom coordinates
			k = 0
			for j in range(num_atoms):
				tmp = data[i].split()
				atoms.append(utils.Atom(tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]) ))
				atoms[-1].fx, atoms[-1].fy, atoms[-1].fz = -grad[k], -grad[k+1], -grad[k+2]
				i += 1
				k += 3
			break
		
	return atoms, energy