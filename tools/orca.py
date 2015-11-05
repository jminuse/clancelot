from merlin import *
from subprocess import Popen

# A function to run an Orca DFT Simulation
def job(run_name, route, atoms=[], extra_section='', grad=False, queue=None, procs=1, charge_and_multiplicity='0 1', title='', blurb=None, force=False, previous=None, neb=[False,None,None,None], mem=None, lint=False):
	# Generate the orca input file
	os.system('mkdir -p orca/%s' % run_name)
	os.chdir('orca/%s' % run_name)

	# If running on system with more than one core, tell orca
	if procs > 1:
		#route += ' PAL%d' % procs
		extra_section = '%pal nprocs '+str(procs)+' end\n' + extra_section.strip()

	# If desiring .orca.engrad output, tell orca
	if grad:
		extra_section = extra_section.strip() + '''\n%method
 RunTyp Gradient
 end'''

 	# One can specify how much memory they want (in MB) per core
 	if mem is not None:
 		extra_section = extra_section.strip() + '\n%maxcore ' + str(mem).strip()

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
	elif queue=='debug':
		print 'Would run', run_name, charge_and_multiplicity
	else:
		NBS = '''#!/bin/bash
##NBS-name: "%s"
##NBS-nproc: %d
##NBS-queue: "%s"

export PATH=/fs/europa/g_pc/ompi_1_6_5/bin:/fs/europa/g_pc/orca_3_0_3_linux_x86-64:$PATH
export LD_LIBRARY_PATH=/fs/europa/g_pc/ompi_1_6_5/lib:$LD_LIBRARY_PATH

/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca %s.orca > %s.out 2>&1
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
def parse_atoms(input_file, get_atoms=True, get_energy=True, get_charges=False, charge_type='MULLIKEN', get_time=False, get_bandgap=False, check_convergence=False, check_converged=False, parse_all=False):
	if not os.path.isfile('orca/%s/%s.out' % (input_file,input_file)):
		print("Error - No output file exists for orca sim %s." % input_file)
		sys.exit()
	
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

		if hold.rfind(s) != -1:
			hold = hold[hold.rfind(s):]
			b = hold[:hold.find('\n\n')].split('\n')[2:-1]
			for a in b:
				a = a.split()
				charges.append([a[1].split(':')[0],float(a[-1])])
		else:
			charges = None
	
	# Get Total Simulation Time
	if get_time:
		hold, times = data, []
		s = 'TOTAL RUN TIME'
		if hold.rfind(s) != -1:
			hold = hold[hold.rfind(s):]
			hold = hold[:hold.find('\n')].split()
			time = float(hold[3])*3600*24 + float(hold[5])*3600 + float(hold[7])*60 + float(hold[9]) + float(hold[11])/1000.0
			times.append(time)
		else:
			times.append(float('NaN'))

	if get_bandgap:
		hold, bandgap = data, []
		s = 'ORBITAL ENERGIES'
		while hold.find(s) != -1:
			hold = hold[hold.find(s) + len(s):]
			tmp = hold[:hold.replace('\n\n','\t\t',1).find('\n\n')].split('\n')[4:]
			for i,t in enumerate(tmp):
				t = t.split()
				if float(t[1]) == 0:
					if i == 0:
						print("Error in calculating bandgap. Lowest energy orbital is empty.")
						sys.exit()
					bandgap.append(float(t[2]) - float(tp[2]))
					break
				tp = t
			hold = hold[hold.find('\n'):]
		if not parse_all: bandgap = bandgap[-1]

	if check_convergence:
		hold, convergence = data, []
		s = 'Geometry convergence'
		if hold.rfind(s) != -1:
			hold = hold[hold.find(s)+len(s):]
			tmp = hold[:hold.find('Max(Bonds)')].split('\n')[3:-2]
			tmp_convergence = []
			for a in tmp:
				a = a.split()
				tmp_convergence.append([' '.join(a[:2]),float(a[2]),float(a[3]), a[4]])
			convergence.append(tmp_convergence)
		else:
			convergence = None

	if check_converged:
		hold, converged = data, False
		s = 'SCF CONVERGED AFTER'
		if hold.find(s) != -1:
			converged = True

	results = []
	if get_atoms: results.append(atoms)
	if get_energy: results.append(energies)
	if get_charges: results.append(charges)
	if get_time: results.append(times)
	if get_bandgap: results.append(bandgap)
	if check_convergence: results.append(convergence)
	if check_converged: results.append(converged)

	while type(results) == list and len(results) == 1:
		results = results[0]
	return results

# Simplified calls to parse_atoms
def atoms(input_file, parse_all=True):
	return parse_atoms(input_file, get_atoms=True, get_energy=False, get_charges=False, get_time=False, get_bandgap=False, check_convergence=False, check_converged=False, parse_all=parse_all)
def energies(input_file, parse_all=True):
	return parse_atoms(input_file, get_atoms=False, get_energy=True, get_charges=False, get_time=False, get_bandgap=False, check_convergence=False, check_converged=False, parse_all=parse_all)
def charges(input_file, parse_all=True):
	return parse_atoms(input_file, get_atoms=False, get_energy=False, get_charges=True, get_time=False, get_bandgap=False, check_convergence=False, check_converged=False, parse_all=parse_all)
def times(input_file, parse_all=True):
	return parse_atoms(input_file, get_atoms=False, get_energy=False, get_charges=False, get_time=True, get_bandgap=False, check_convergence=False, check_converged=False, parse_all=parse_all)
def bandgap(input_file, parse_all=True):
	return parse_atoms(input_file, get_atoms=False, get_energy=False, get_charges=False, get_time=False, get_bandgap=True, check_convergence=False, check_converged=False, parse_all=parse_all)
def convergence(input_file, parse_all=True):
	return parse_atoms(input_file, get_atoms=False, get_energy=False, get_charges=False, get_time=False, get_bandgap=False, check_convergence=True, check_converged=False, parse_all=parse_all)
def converged(input_file, parse_all=True):
	return parse_atoms(input_file, get_atoms=False, get_energy=False, get_charges=False, get_time=False, get_bandgap=False, check_convergence=False, check_converged=True, parse_all=parse_all)

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

def read(input_file):
	data = utils.DFT_out(input_file, 'orca')

	data.frames = atoms(input_file, parse_all=True)
	if type(data.frames) == list and type(data.frames[0]) != list: data.frames = [data.frames]
	data.atoms = data.frames[-1] if type(data.frames)==list and type(data.frames[0])==list else data.frames
	data.energies = energies(input_file, parse_all=True)
	data.charges_MULLIKEN = parse_atoms(input_file, get_atoms=False, get_energy=False, get_charges=True, charge_type='MULLIKEN', get_time=False, get_bandgap=False, check_convergence=False, parse_all=True)
	data.charges_LOEWDIN = parse_atoms(input_file, get_atoms=False, get_energy=False, get_charges=True, charge_type='LOEWDIN', get_time=False, get_bandgap=False, check_convergence=False, parse_all=True)
	data.charges = data.charges_MULLIKEN if data.charges_MULLIKEN is not None else data.charges_LOEWDIN
	data.convergence = convergence(input_file, parse_all=True)
	data.converged = converged(input_file, parse_all=True)
	data.time = times(input_file, parse_all=True)
	data.bandgap = bandgap(input_file, parse_all=True)


	return data
