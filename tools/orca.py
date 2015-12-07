from merlin import *
from subprocess import Popen
from copy import deepcopy
from shutil import copyfile

def read(input_file):
	# Check file exists, and open
	if not os.path.isfile('orca/%s/%s.out' % (input_file,input_file)):
		print("Error - No output file exists for orca sim %s." % input_file)
		sys.exit()
	data = open('orca/%s/%s.out' % (input_file,input_file),'r').read()

	# Get all the positions
	hold, frames = data, []
	s = 'CARTESIAN COORDINATES (ANGSTROEM)'
	while hold.find(s) != -1:
		hold = hold[hold.find(s)+len(s):]
		tmp = hold[:hold.find('\n\n')].split('\n')[2:]
		tmp_atoms = []
		for a in tmp:
			a = a.split()
			tmp_atoms.append(utils.Atom(a[0],float(a[1]),float(a[2]),float(a[3])))
		frames.append(tmp_atoms)

	if len(frames) > 0:
		atoms = frames[-1]
	else:
		atoms = None

	# Get all the energies
	hold, energies = data, []
	s = 'FINAL SINGLE POINT ENERGY'
	while hold.find(s) != -1:
		hold = hold[hold.find(s):]
		tmp = hold[:hold.find('\n')].split()[-1]
		energies.append(float(tmp))
		hold = hold[hold.find('\n'):]
	
	if len(energies) > 0:
		energy = energies[-1]
	else:
		energy = None

	# Get charges
	hold, charges_MULLIKEN = data, []
	s = 'MULLIKEN ATOMIC CHARGES'
	if hold.rfind(s) != -1:
		hold = hold[hold.rfind(s):]
		b = hold[:hold.find('\n\n')].split('\n')[2:-1]
		for a in b:
			a = a.split()
			charges_MULLIKEN.append([a[1].split(':')[0],float(a[-1])])
	else:
		charges_MULLIKEN = None

	hold, charges_LOEWDIN = data, []
	s = 'LOEWDIN ATOMIC CHARGES'
	if hold.rfind(s) != -1:
		hold = hold[hold.rfind(s):]
		b = hold[:hold.find('\n\n')].split('\n')[2:]
		for a in b:
			a = a.split()
			charges_LOEWDIN.append([a[1].split(':')[0],float(a[-1])])
	else:
		charges_LOEWDIN = None

	# Get Total Simulation Time
	hold = data
	s = 'TOTAL RUN TIME'
	if hold.rfind(s) != -1:
		hold = hold[hold.rfind(s):]
		hold = hold[:hold.find('\n')].split()
		time = float(hold[3])*3600*24 + float(hold[5])*3600 + float(hold[7])*60 + float(hold[9]) + float(hold[11])/1000.0
	else:
		time = float('NaN')

	hold, bandgaps = data, []
	s = 'ORBITAL ENERGIES'
	while hold.find(s) != -1:
		hold = hold[hold.find(s) + len(s):]
		tmp = hold[:hold.replace('\n\n','\t\t',1).find('\n\n')].split('\n')[4:]
		for i,t in enumerate(tmp):
			t = t.split()
			if float(t[1]) == 0:
				if i == 0:
					print("Error in calculating bandgaps. Lowest energy orbital is empty.")
					sys.exit()
				bandgaps.append(float(t[2]) - float(tp[2]))
				break
			tp = t
		hold = hold[hold.find('\n'):]
	
	if len(bandgaps) > 0:
		bandgap = bandgaps[-1]
	else:
		bandgap = None

	hold, convergence = data, []
	s = 'Geometry convergence'
	if hold.rfind(s) != -1:
		hold = hold[hold.find(s)+len(s):]
		tmp = hold[:hold.find('Max(Bonds)')].split('\n')[3:-2]
		convergence = []
		for a in tmp:
			a = a.split()
			convergence.append([' '.join(a[:2]),float(a[2]),float(a[3]), a[4]])
	else:
		convergence = None

	hold, converged = data, False
	s = 'SCF CONVERGED AFTER'
	if hold.find(s) != -1:
		converged = True

	data = utils.DFT_out(input_file, 'orca')

	data.frames = frames
	data.atoms = atoms
	data.energies = energies
	data.energy = energy
	data.charges_MULLIKEN = charges_MULLIKEN
	data.charges_LOEWDIN = charges_LOEWDIN
	data.charges = deepcopy(charges_MULLIKEN)
	data.convergence = convergence
	data.converged = converged
	data.time = time
	data.bandgaps = bandgaps
	data.bandgap = bandgap

	return data

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

# A function to run an Orca DFT Simulation
def job(run_name, route, atoms=[], extra_section='', grad=False, queue=None, procs=1, charge_and_multiplicity='0 1', previous=None, mem=None):
	# Generate the orca input file
	os.system('mkdir -p orca/%s' % run_name)
	os.chdir('orca/%s' % run_name)

	# If running on system with more than one core, tell orca
	if procs > 1:
		#route += ' PAL%d' % procs
		extra_section = '%pal nprocs '+str(procs)+' end\n' + extra_section.strip()

	# If no atoms were specified, but previous was, grab lowest energy atoms
	if previous is not None and atoms == []:
		pwd = os.getcwd()
		os.chdir('../../')
		hold = read(previous)
		os.chdir(pwd)
		index = hold.energies.index(min(hold.energies))
		atoms = hold.frames[index]
		hold = None

	# If desiring .orca.engrad output, tell orca
	if grad:
		extra_section = extra_section.strip() + '''\n%method
 RunTyp Gradient
 end'''

 	# One can specify how much memory they want (in MB) per core
 	if mem is not None:
 		extra_section = extra_section.strip() + '\n%maxcore ' + str(mem).strip()

 	# Add moread if previous
 	if previous is not None:
 		route = route.strip() + ' MORead'
 		PATH = '../'+previous+'/'+previous+'.orca.gbw'
 		if not os.path.isfile(PATH):
 			print("Error - Previous run %s does not have a .gbw file." % previous)
 			sys.exit()
 		extra_section = extra_section.strip() + '\n%moinp "' + PATH + '"' 

 	# ---------------------------------------------------------------------------
 	# NO MORE CHANGES TO EXTRA_SECTION AFTER THIS!-------------------------------
 	# ---------------------------------------------------------------------------

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

	# Copy run script
	fname = sys.argv[0]
	if '/' in fname: fname = fname.split('/')[-1]
	try:
		copyfile('../../%s' % fname,fname)
	except IOError:
		# Submitted a job oddly enough that sys.argv[0] is not the original python file name, so don't do this
		pass

	# Return to the appropriate directory
	os.chdir('../..')

	if queue is None:
		return process_handle
	else:
		return utils.Job(run_name)

