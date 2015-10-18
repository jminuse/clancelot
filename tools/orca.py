from merlin import *

# A function to run an Orca DFT Simulation
def job(run_name, route, atoms=[], extra_section='', queue=None, procs=1, charge_and_multiplicity='0 1', title='', blurb=None, force=False, previous=None, neb=[False,None,None,None], lint=False):
	# Generate the orca input file
	os.system('mkdir -p orca/%s' % run_name)
	os.chdir('orca/%s' % run_name)

	# If running on system with more than one core, tell orca
	if procs > 1:
		route += ' PAL%d' % procs

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
		os.system('/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca %s.orca > %s.out &' % (run_name,run_name))
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

# A function to parse and Orca DFT output file (assumes by default a .out file format)
def parse_atoms(input_file, get_atoms=True, get_energy=True, get_charges=False, charge_type='MULLIKEN', get_time=False, check_convergence=False, parse_all=False):
	# This function returns data as:
	#	atoms, energy, charges, time, convergence
	# If parse_all is on, everything is in lists containing all information from the simulation
	data = open('orca/%s/%s.out' % (input_file,input_file),'r').read()
	atoms, energies, charges, times, convergence = [], [], [], [], []
	tmp_atoms, tmp_energy, tmp_time, tmp_convergence, chk_done, incomplete = [], None, None, [], False, False
	while data.find('CARTESIAN COORDINATES (ANGSTROEM)') != -1:
		tmp_atoms, tmp_energy, tmp_times, tmp_convergence = [], None, None, []
		
		s = 'CARTESIAN COORDINATES (ANGSTROEM)'
		if get_atoms and data.find(s)!=-1:
			data = data[data.find(s)+len(s):]
			b = data[:data.find('\n\n')].split('\n')[2:]
			for a in b:
				a = a.split()
				tmp_atoms.append(utils.Atom(a[0],float(a[1]),float(a[2]),float(a[3])))
			tmp_atoms_hold = tmp_atoms
		elif get_atoms:
			tmp_atoms = tmp_atoms_hold
			incomplete = True

		s = 'Total SCF time'
		if get_time and data.find(s)!=-1:
			data = data[data.find(s):]
			b = data[:data.find('\n')].split()
			data = data[data.find('\n'):]
			tmp_times = float(b[3])*3600*24 + float(b[5])*3600 + float(b[7])*60 + float(b[9])
			tmp_times_hold = tmp_times
		elif get_time:
			tmp_times = tmp_times_hold
			incomplete = True
		
		s = 'FINAL SINGLE POINT ENERGY'
		if get_energy and data.find(s)!=-1:
			data = data[data.find(s):]
			b = data[:data.find('\n')].split()[-1]
			data = data[data.find('\n'):]
			tmp_energy = float(b)
			tmp_energy_hold = tmp_energy
		elif get_energy:
			tmp_energy = tmp_energy_hold
			incomplete = True
		
		if check_convergence and data.find('Geometry convergence') != -1:
			s = 'Geometry convergence'
			data = data[data.find(s)+len(s):]
			b = data[:data.find('Max(Bonds)')].split('\n')[3:-2]
			for a in b:
				a = a.split()
				tmp_convergence.append([' '.join(a[:2]),float(a[2]),float(a[3]), a[4]=='YES'])
			tmp_convergence_hold = tmp_convergence
		elif check_convergence and not chk_done:
			tmp_convergence = tmp_convergence_hold
			chk_done = True
		elif chk_done:
			print("Error - Could not get the geometry convergence value")
			sys.exit()

		if parse_all:
			if get_atoms: atoms.append(tmp_atoms)
			if get_time: times.append(tmp_time)
			if get_energy: energies.append(tmp_energy)
			if check_convergence: convergence.append(tmp_convergence)
		if incomplete:
			break

	if get_charges:
		if charge_type.strip().upper() not in ['MULLIKEN','LOEWDIN']:
			print("Error - Requested %s charge, but not a valid charge to request.  Choose from MULLIKEN or LOEWDIN." % charge_type.strip().upper())
			sys.exit()
		s = charge_type.strip().upper()+' ATOMIC CHARGES'
		data = open('orca/%s/%s.out' % (input_file,input_file),'r').read()
		data = data[data.rfind(s):]
		b = data[:data.find('\n\n')].split('\n')[2:-1]
		for a in b:
			a = a.split()
			charges.append([a[1],float(a[3])])
	
	results = []
	if get_atoms: results.append(atoms)
	if get_energy: results.append(energies)
	if get_charges: results.append(charges)
	if get_time: results.append(times)
	if check_convergence: results.append(convergence)
	return results
