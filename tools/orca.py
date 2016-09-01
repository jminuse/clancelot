from merlin import *
from subprocess import Popen
from copy import deepcopy
from shutil import copyfile
from time import sleep
import re

def read(input_file):
	# Check file exists, and open
	if input_file.startswith('/'): #allow absolute paths as filenames
		input_path = input_file
	elif os.path.isfile(input_file):
		input_path = input_file
	else:
		input_path = 'orca/%s/%s.out' % (input_file,input_file)
	if not os.path.isfile(input_path):
		raise IOError('Expected orca output file does not exist at %s' % (input_path))
		sys.exit()
	data = open(input_path,'r').read()
	data_lines = data.splitlines()

	# Get the route line
	try:
		route = [line[5:] for line in data_lines if line.startswith('|  1>')][0]
	except IndexError:
		raise IOError('Could not find route line in %s: job most likely crashed.' % input_path)

	# Get all the energies
	energies = [float(e) for e in re.findall('FINAL SINGLE POINT ENERGY +(\S+)', data)]
	
	if len(energies) > 0:
		energy = min(energies)
	else:
		energy = None
		
	# Get all the positions
	section, frames = data, []
	s = 'CARTESIAN COORDINATES (ANGSTROEM)'
	while s in section:
		section = section[section.find(s)+len(s):]
		atom_block = section[:section.find('\n\n')].split('\n')[2:]
		frame = []
		for line in atom_block:
			a = line.split()
			frame.append(utils.Atom(a[0],float(a[1]),float(a[2]),float(a[3])))
		frames.append(frame)

	if frames:
		if energy is None:
			atoms = frames[-1]
		else:
			atoms = frames[ energies.index(energy) ]
	else:
		atoms = None

	# Get all the gradients if CARTESIAN GRADIENTS is in the file.
	# Else, if MP2 gradients is in the file, grab the last gradient
	s_gradient = "CARTESIAN GRADIENT"
	s_gradient_2 = "The final MP2 gradient"
	section, gradients = data, []
	if s_gradient in section:
		s = s_gradient
	elif s_gradient_2 in section:
		s = s_gradient_2
	else:
		s, gradients = None, None
	if s is not None:
		while s in section:
			gradient = []
			if s == s_gradient:
				grad_block = section[section.find(s_gradient):].split("\n\n")[1].split("\n")
				gradient = []
				for line in grad_block:
					a = line.split()
					gradient.append([float(b) for b in a[3:]])
			elif s == s_gradient_2:
				grad_block = section[section.find(s_gradient_2):].split("\n\n")[0].split("\n")[1:]
				gradient = []
				for line in grad_block:
					a = line.split()
					gradient.append([float(b) for b in a[1:]])
			section = section[section.find(s)+len(s):]
			gradients.append(gradient)

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
		for a, charge in zip(atoms, charges_LOEWDIN):
			a.charge = charge[1]
	else:
		charges_LOEWDIN = None

	hold, charges_CHELPG = data, []
	s = 'CHELPG Charges'
	if hold.rfind(s) != -1:
		hold = hold[hold.rfind(s):]
		b = hold[:hold.find('\n--------------------------------\nTotal charge:')].split('\n')[2:]
		for a in b:
			a = a.split()
			charges_CHELPG.append([a[1].split(':')[0],float(a[-1])])
		for a, charge in zip(atoms, charges_CHELPG):
			a.charge = charge[1]
	else:
		charges_CHELPG = None

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
		hold = hold[hold.rfind(s)+len(s):]

		# Cartesian optimization does not compute Max(Bonds). Instead use a more general '\n\n' if 'Max(Bonds)' cannot be found
		if hold.rfind('Max(Bonds)') != -1:
			tmp = hold[:hold.rfind('Max(Bonds)')].split('\n')[3:-2]
		else:
			tmp = hold[:hold.find('\n\n')].split('\n')[3:-1]

		convergence = []
		for a in tmp:
			a = a.split()
			convergence.append([' '.join(a[:2]),float(a[2]),float(a[3]), a[4]])
	else:
		convergence = None

	hold, converged = data, False
	s1, s2 = 'SCF CONVERGED AFTER', 'OPTIMIZATION RUN DONE'
	if 'opt' in route.lower(): s = s2
	else: s = s1
	if hold.find(s) != -1:
		converged = True

	finished = 'ORCA TERMINATED NORMALLY' in data
		
	warnings = [line for line in data_lines if line.startswith('Warning: ')]

	data = utils.DFT_out(input_file, 'orca')

	data.route = route
	data.frames = frames
	data.atoms = atoms
	data.gradients = gradients
	data.energies = energies
	data.energy = energy
	data.charges_MULLIKEN = charges_MULLIKEN
	data.charges_LOEWDIN = charges_LOEWDIN
	data.charges_CHELPG = charges_CHELPG
	data.charges = deepcopy(charges_MULLIKEN)
	data.convergence = convergence
	data.converged = converged
	data.time = time
	data.bandgaps = bandgaps
	data.bandgap = bandgap
	data.finished = finished
	data.warnings = warnings

	return data

# A function to parse orca.engrad files
def engrad_read(input_file, force='Ha/Bohr', pos='Bohr'):
	if not input_file.endswith('.engrad'):
		input_file = 'orca/%s/%s.orca.engrad' % (input_file,input_file)
	if not os.path.isfile(input_file):
		raise IOError("No engrad file %s exists in %s. Please run simulation with grad=True." % (input_file,os.getcwd()))
	data = open(input_file,'r').read().split('\n')
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
				atoms.append(utils.Atom(tmp[0], units.convert_dist('Bohr',pos,float(tmp[1])),
								units.convert_dist('Bohr',pos,float(tmp[2])),
								units.convert_dist('Bohr',pos,float(tmp[3]))))
				atoms[-1].fx = units.convert('Ha/Bohr',force,-grad[k])
				atoms[-1].fy = units.convert('Ha/Bohr',force,-grad[k+1])
				atoms[-1].fz = units.convert('Ha/Bohr',force,-grad[k+2])
				i += 1
				k += 3
			break
		
	return atoms, energy

# A function to run an Orca DFT Simulation
def job(run_name, route, atoms=[], extra_section='', grad=False, queue=None, procs=1, charge=None, multiplicity=None, charge_and_multiplicity='0 1', previous=None, mem=2000, priority=100, xhost=None):
	if len(run_name) > 31 and queue is not None:
		raise Exception("Job name too long (%d) for NBS. Max character length is 31." % len(run_name))

	# Generate the orca input file
	os.system('mkdir -p orca/%s' % run_name)
	os.chdir('orca/%s' % run_name)

	if '!' not in route:
		route = '! ' + route #orca requires route line to start with "!"

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

 	if xhost is not None:
 		if type(xhost) is str:
 			xhosts = "##NBS-xhost: \"%s\"" % xhost
 		else:
 			xhosts = "##NBS-xhost: " + ", ".join( map(lambda x: '"' + x + '"', xhost) )
 	else:
 		xhosts = ""

 	# Add moread if a previous orca job was provided
 	if previous is not None:
 		route = route.strip() + ' MORead'
 		current_dir = os.getcwd()
 		if previous.startswith('/'): #accept absolute paths
 			previous_path = previous
 		else:
	 		previous_path = current_dir+'/../'+previous+'/'+previous+'.orca.gbw'
	 		if not os.path.isfile(previous_path):
	 			if os.path.isfile(current_dir+'/../'+previous+'/'+previous+'.orca.proc0.gbw'):
	 				previous_path = current_dir+'/../'+previous+'/'+previous+'.orca.proc0.gbw'
 		if not os.path.isfile(previous_path):
 			raise Exception("Previous run does not have a .gbw file at %s." % (previous_path))
 		copyfile(previous_path, 'previous.gbw')
 		extra_section = extra_section.strip() + '\n%moinp "previous.gbw"'
 		
 		
 		# If no atoms were specified, get the atoms from the previous job
		if atoms == []:
			old_results = read(previous_path.replace('.orca.gbw','.out').replace('.orca.proc0.gbw','.out'))
			if " opt " in route.lower() or " sp " in route.lower(): # then take lowest-energy atoms from previous job
				index = old_results.energies.index(min(old_results.energies)) 
				atoms = old_results.frames[index]
			else: #otherwise just take the last atoms
				atoms = old_results.frames[-1]

 	# ---------------------------------------------------------------------------
 	# NO MORE CHANGES TO EXTRA_SECTION AFTER THIS!-------------------------------
 	# ---------------------------------------------------------------------------

	if charge!=None or multiplicity!=None: #if either charge or multiplicity are specified, use that
		if charge==None: charge = 0 #default neutral charge
		if multiplicity==None: multiplicity = 1 #default singlet state
		charge_and_multiplicity = '%d %d' % (charge, multiplicity)
	# Get input for orca formatted correctly
	inp = route.strip()+'\n'+extra_section.strip()+'\n*xyz '+charge_and_multiplicity+'\n'
	for a in atoms:
		inp += '%s %f %f %f %s\n' % (a.element, a.x, a.y, a.z, (a.extra if hasattr(a,'extra') else ''))
	inp += '*\n'

	# Write the orca file
	f = open(run_name+'.orca', 'w')
	f.write(inp)
	f.close()

	# Run the simulation
	if queue is None:
		process_handle = Popen('/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca %s.orca > %s.out' % (run_name, run_name), shell=True)
	elif queue=='debug':
		print 'Would run', run_name
	else: #details copied from g09sub
		
		NBS = '''#!/bin/sh
##NBS-fdisk: 8192
##NBS-fmemory: '''+str(mem)+'''
##NBS-fswap: 8192
##NBS-unique: yes
##NBS-sandbox: yes
##NBS-tmp-sandbox: yes
##NBS-name: '''+run_name+'''
##NBS-nproc: '''+str(procs)+'''
##NBS-queue: '''+queue+'''
##NBS-priority: '''+str(priority)+'''
'''+xhosts+'''

##NBS-input: *.orca
'''+('##NBS-input: previous.gbw' if os.path.exists('previous.gbw') else '')+'''

##NBS-output: *.xyz -overwrite
##NBS-output: *.trj -overwrite
##NBS-output: *.gbw -overwrite
##NBS-output: *.engrad -overwrite
##NBS-output: *.prop -overwrite
##NBS-output: *.opt -overwrite

export PATH=/fs/europa/g_pc/ompi_1_6_5/bin:/fs/europa/g_pc/orca_3_0_3_linux_x86-64:$PATH
export LD_LIBRARY_PATH=/fs/europa/g_pc/ompi_1_6_5/lib:$LD_LIBRARY_PATH

/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca '''+run_name+'''.orca > '''+(os.getcwd()+'/'+run_name)+'''.out

cp /tmp/*/'''+run_name+'''.orca.engrad '''+(os.getcwd()+'/')+'''
cp /tmp/*/'''+run_name+'''.orca.gbw '''+(os.getcwd()+'/')+'''

touch '''+run_name+'''.orca.xyz
touch '''+run_name+'''.orca.trj
touch '''+run_name+'''.orca.gbw
touch '''+run_name+'''.orca.engrad
touch '''+run_name+'''.orca.prop
touch '''+run_name+'''.orca.opt

''' # The use of "touch" above is a patch for an NBS bug, where files cannot be missing or they prevent the others from being copied. 
		f = open(run_name+'.nbs', 'w')
		f.write(NBS)
		f.close()
		os.system('jsub %s.nbs -prop orca' % run_name)
		sleep(0.5)

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

