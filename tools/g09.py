import os, string, sys, re, shutil, copy
from subprocess import Popen
import utils, log, files

def job(run_name, route, atoms=[], extra_section='', queue='batch', procs=1, charge_and_multiplicity='0,1', title='run by gaussian.py', blurb=None, watch=False, eRec=True, force=False, previous=None):
	log.chk_gaussian(run_name,force=force) # Checks if run exists
	head = '#N '+route+'\n\n'+title+'\n\n'+charge_and_multiplicity+'\n' # Header for inp file

	# Setup list of atoms for inp file --------------------------------------------------------------------------------------------
	if atoms and type(atoms[0])==type([]): #multiple lists of atoms (e.g. transistion state calculation)
		xyz = (title+'\n\n0,1\n').join([('\n'.join( [( "%s %f %f %f" % (a.element, a.x, a.y, a.z) ) for a in atom_list] ) + '\n\n') for atom_list in atoms])
	else: #single list of atoms
		if 'oniom' in route.lower():
			xyz = '\n'.join( [( "%s 0 %f %f %f %s" % (a.element, a.x, a.y, a.z, a.layer) ) for a in atoms] ) + '\n\n'
		elif 'counterpoise' in route.lower():
			xyz = '\n'.join( [( "%s(Fragment=%d) %f %f %f" % (a.element, a.fragment, a.x, a.y, a.z) ) for a in atoms] ) + '\n\n'
		elif atoms:
			xyz = '\n'.join( [( "%s %f %f %f" % (a.element, a.x, a.y, a.z) ) for a in atoms] ) + '\n\n'
		else:
			xyz = '\n'

	# cd into gaussian directory, set up .inp file, run simulation, cd out -------------------------------------------------------------------------
	os.chdir('gaussian')
	if queue is not None: #run on queue
		with open(run_name+'.inp', 'w') as inp:
			inp.write(head+xyz+extra_section)
		if previous:	
			shutil.copyfile(previous+'.chk', run_name+'.chk')
		os.system('g09sub '+run_name+' -chk -queue '+queue+((' -nproc '+str(procs)+' ') if procs else '')+' ') #-xhost sys_eei sys_icse
	else: #run not on queue; will hang process until complete
		with open(run_name+'.inp', 'w') as inp:
			csh = '''setenv g09root /usr/local/gaussian/g09d01
source $g09root/g09/bsd/g09.login
g09 <<END > '''+run_name+'''.log
%NProcShared=1
%RWF=/tmp/
%Chk='''+run_name+'''.chk
%Mem=25MW
'''
			inp.write(csh+head+xyz+extra_section+'\neof\nrm /tmp/*.rwf')
		if previous:	
			shutil.copyfile(previous+'.chk', run_name+'.chk')
		process_handle = Popen('/bin/csh %s.inp' % run_name, shell=True)
	shutil.copyfile('../'+sys.argv[0], run_name+'.py')
	os.chdir('..')

	log.put_gaussian(run_name,route,extra_section,blurb,eRec,force) # Places info in log if simulation is run
	if queue is None:
		return process_handle
	else:
		return utils.Job(run_name)
	# ----------------------------------------------------------------------------------------------------------------------------

def restart_job(old_run_name, job_type='ChkBasis Opt=Restart', queue='batch', procs=None):
	run_name = old_run_name+'r'
	os.chdir('gaussian')
	shutil.copyfile(old_run_name+'.chk', run_name+'.chk')
	with open(run_name+'.inp', 'w') as inp:
		inp.write('#t '+job_type+'\n\nrun by gaussian.py\n\n')
	os.system('g09sub '+run_name+' -chk -queue '+queue+((' -nproc '+str(procs)+' ') if procs else '')+' -xhost sys_eei sys_icse')
	os.chdir('..')

def parse_atoms(input_file, get_atoms=True, get_energy=True, check_convergence=True, get_time=False, counterpoise=False):
	if input_file[-4:] != '.log':
		input_file = 'gaussian/'+input_file+'.log'
	
	contents = open(input_file).read()
	if check_convergence and get_energy and 'Normal termination of Gaussian 09' not in contents:
		return None

	if 'Summary of Optimized Potential Surface Scan' in contents:
		end_section = contents[contents.rindex('Summary of Optimized Potential Surface Scan'):]
		energy_lines = re.findall('Eigenvalues -- ([^\\n]+)', end_section)
		energy = [float(s) for line in energy_lines for s in re.findall('-[\d]+\.[\d]+', line)]
		
		minima = re.split('Stationary point found', contents)
		atoms = []
		for m in minima[1:]:
			coordinates = m.index('Coordinates (Angstroms)')
			
			start = m.index('---\n', coordinates)+4
			end = m.index('\n ---', start)
			atoms.append([])
			for line in m[start:end].splitlines():
				columns = line.split()
				element = columns[1]
				x,y,z = [float(s) for s in columns[3:6]]
				atoms[-1].append( utils.Atom(element=utils.elements_by_atomic_number[int(columns[1])], x=x, y=y, z=z, index=len(atoms[-1])+1) )
		
		if get_energy:
			return energy, atoms
		
	elif get_energy:
		if ' MP2/' in contents: # MP2 files don't have just SCF energy
			energy = float(re.findall('EUMP2 = +(\S+)', contents)[-1].replace('D','e'))
		elif ' CCSD/' in contents:
			energy = float(re.findall('E\(CORR\)= +(\S+)', contents)[-1])
		else:
			if not counterpoise:
				try:
					energy_line = contents[contents.rindex('SCF Done'):contents.index('\n', contents.rindex('SCF Done'))]
				except ValueError:
					raise Exception('No SCF for '+input_file)
				energy = float(re.search('SCF Done: +\S+ += +(\S+)', energy_line).group(1))
			else:
				energy = float(re.findall('Counterpoise: corrected energy = +(\S+)', contents)[-1])
	
	if get_time:
		m = re.search('Job cpu time: +(\S+) +days +(\S+) +hours +(\S+) +minutes +(\S+) +seconds', contents)
		time = float(m.group(1))*24*60*60 + float(m.group(2))*60*60 + float(m.group(3))*60 + float(m.group(4))
	
	if get_energy and not get_atoms:
		if get_time:
			return energy, time
		else:
			return energy

	#get coordinates
	try:
		last_coordinates = contents.rindex('Input orientation:')
		last_coordinates = contents.index('Coordinates (Angstroms)', last_coordinates)
	except ValueError:
		last_coordinates = contents.index('Coordinates (Angstroms)')
	start = contents.index('---\n', last_coordinates)+4
	end = contents.index('\n ---', start)
	atoms = []
	for line in contents[start:end].splitlines():
		columns = line.split()
		element = columns[1]
		x,y,z = [float(s) for s in columns[3:6]]
		atoms.append( utils.Atom(element=utils.elements_by_atomic_number[int(columns[1])], x=x, y=y, z=z, index=len(atoms)+1) )
	
	#get forces
	if 'Forces (Hartrees/Bohr)' in contents:
		last_forces = contents.rindex('Forces (Hartrees/Bohr)')
		start = contents.index('---\n', last_forces)+4
		end = contents.index('\n ---', start)
		for i,line in enumerate(contents[start:end].splitlines()):
			columns = line.split()
			atoms[i].fx, atoms[i].fy, atoms[i].fz = [float(s) for s in columns[2:5]]
	
	#return the appropriate values
	if get_time:
		if get_atoms:
			return energy, atoms, time
		else:
			return energy, time
	if get_energy:
		return energy, atoms
	else:
		return atoms
		
def atoms(input_file, check=False):
	return parse_atoms(input_file, get_atoms=True, get_energy=False, check_convergence=check, get_time=False, counterpoise=False)

def parse_all(input_file):
	contents = open(input_file).read()
	time = None
	if 'Normal termination of Gaussian 09' not in contents:
		pass
	else:
		m = re.search('Job cpu time: +(\S+) +days +(\S+) +hours +(\S+) +minutes +(\S+) +seconds', contents)
		time = float(m.group(1))*24*60*60 + float(m.group(2))*60*60 + float(m.group(3))*60 + float(m.group(4))

	energies = []
	atom_frames = []
	start = 0
	while True:
		try: #match energy
			start = contents.index('SCF Done', start)
			energy_this_step = float( re.search('SCF Done: +\S+ += +(\S+)', contents[start:]).group(1) )
			input_orientation = contents.find('Input orientation:', start)
			if input_orientation >= 0:
				start = input_orientation
			next_coordinates = contents.index('Coordinates (Angstroms)', start)
		except: break
		start = contents.index('---\n', next_coordinates)+4
		end = contents.index('\n ---', start)
		lines = contents[start:end].splitlines()
		start = end

		atoms = []
		for line in lines:
			columns = line.split()
			element = columns[1]
			x,y,z = columns[3:6]
			atoms.append( utils.Atom(element=element, x=float(x), y=float(y), z=float(z)) )
		atom_frames.append(atoms)
		energies.append(energy_this_step)

	return energies, atom_frames, time
	
def parse_scan(input_file):
	contents = open(input_file).read()
	if 'Normal termination of Gaussian 09' not in contents:
		return None
	scan_steps = contents.split('on scan point')
	energy_list = []
	atoms_list = []
	
	scan_steps = [ scan_steps[i] for i in range(1,len(scan_steps)-1) if scan_steps[i][:10].split()[0]!=scan_steps[i+1][:10].split()[0] ]
	#print [int(s[:10].split()[0]) for s in scan_steps]
	#print len(scan_steps)
	
	for scan_step in scan_steps:
		energy_line = scan_step[scan_step.rindex('SCF Done'):scan_step.index('\n', scan_step.rindex('SCF Done'))]
		energy = float(re.search('SCF Done: +\S+ += +(\S+)', energy_line).group(1))

		last_coordinates = scan_step.rindex('Coordinates (Angstroms)')

		start = scan_step.index('---\n', last_coordinates)+4
		end = scan_step.index('\n ---', start)
		atoms = []
		for line in scan_step[start:end].splitlines():
			columns = line.split()
			element = columns[1]
			x,y,z = [float(s) for s in columns[3:6]]
			atoms.append( utils.Atom(element=utils.elements_by_atomic_number[int(columns[1])], x=x, y=y, z=z) )
		energy_list.append(energy)
		atoms_list.append(atoms)
	return energy_list, atoms_list
	
def parse_chelpg(input_file):
	with open(input_file) as inp:
		contents = inp.read()
	if 'Normal termination of Gaussian 09' not in contents:
		return None
	
	start = contents.rindex('Fitting point charges to electrostatic potential')
	end = contents.index('-----------------', start)
	charges = []
	for line in contents[start:end].splitlines():
		columns = line.split()
		if len(columns)==3:
			charges.append( float(columns[2]) )
	return charges

def neb(name, states, theory, extra_section='', queue=None, spring_atoms=None, k=0.1837, fit_rigid=True): #Nudged Elastic Band. k for VASP is 5 eV/Angstrom, ie 0.1837 Hartree/Angstrom. 
#Cite NEB: http://scitation.aip.org/content/aip/journal/jcp/113/22/10.1063/1.1323224
	from scipy.optimize import minimize
	import numpy as np
	#set which atoms will be affected by virtual springs
	if not spring_atoms:#if not given, select all
		spring_atoms = range(len(states[0]))
	elif type(spring_atoms)==str: #a list of element names
		elements = spring_atoms.split()
		spring_atoms = [i for i,a in enumerate(states[0]) if a.element in elements]
	#class to contain working variables
	class NEB:
		name, states, theory, k = None, None, None, None
		error, forces = None, None
		step = 0
		def __init__(self, name, states, theory, k=1e-2):
			NEB.name = name
			NEB.states = states
			NEB.theory = theory
			NEB.k = k
			
			if fit_rigid:
				#center all states around spring-held atoms
				for s in states:
					center_x = sum([a.x for i,a in enumerate(s) if i in spring_atoms])/len(spring_atoms)
					center_y = sum([a.y for i,a in enumerate(s) if i in spring_atoms])/len(spring_atoms)
					center_z = sum([a.z for i,a in enumerate(s) if i in spring_atoms])/len(spring_atoms)
					for a in s:
						a.x -= center_x
						a.y -= center_y
						a.z -= center_z
				#rotate all states to be as similar to their neighbors as possible
				from scipy.linalg import orthogonal_procrustes
				for i in range(1,len(states)): #rotate all states to optimal alignment
					#only count spring-held atoms for finding alignment
					spring_atoms_1 = [(a.x,a.y,a.z) for j,a in enumerate(states[i]) if j in spring_atoms]
					spring_atoms_2 = [(a.x,a.y,a.z) for j,a in enumerate(states[i-1]) if j in spring_atoms]
					rotation = orthogonal_procrustes(spring_atoms_1,spring_atoms_2)[0]
					#rotate all atoms into alignment
					for a in states[i]:
						a.x,a.y,a.z = utils.matvec(rotation, (a.x,a.y,a.z))
	
			#load initial coordinates into flat array for optimizer
			NEB.coords_start = []
			for s in states[1:-1]:
				for a in s:
					NEB.coords_start += [a.x, a.y, a.z]
	
		@staticmethod
		def calculate(coords):
			coord_count = 0
			for s in NEB.states[1:-1]:
				for a in s:
					a.x, a.y, a.z = coords[coord_count], coords[coord_count+1], coords[coord_count+2]
					coord_count += 3
			#start DFT jobs
			running_jobs = []
			for i,state in enumerate(NEB.states):
				if NEB.step>0:
					if (i==0 or i==len(NEB.states)-1): #only perform evaluations of endpoints on first step
						continue
					guess = ' Guess=Read'
				else: guess = '' #no previous guess for first step
				running_jobs.append( job('%s-%d-%d'%(NEB.name,NEB.step,i), NEB.theory+' Force'+guess, state, queue=queue, force=True, previous=('%s-%d-%d'%(NEB.name,NEB.step-1,i)) if NEB.step>0 else None, extra_section=extra_section) )
			#wait for jobs to finish
			for j in running_jobs: j.wait()
			#get forces and energies from DFT calculations
			energies = []
			for i,state in enumerate(NEB.states):
				try:
					if (i==0 or i==len(NEB.states)-1): # if at an endpoint, just use first step's result
						step_to_use = 0
					else:
						step_to_use = NEB.step
					new_energy, new_atoms = parse_atoms('%s-%d-%d' % (NEB.name, step_to_use, i), check_convergence=False)
				except:
					print "Unexpected error in 'parse_atoms':", sys.exc_info()[0]
					print 'Job failed: %s-%d-%d'%(NEB.name,NEB.step,i); exit()
				energies.append(new_energy)
				for a,b in zip(state, new_atoms):
					a.fx = b.fx; a.fy = b.fy; a.fz = b.fz
			V = copy.deepcopy(energies) # V = potential energy from DFT. energies = V+springs
			#add spring forces to atoms
			for i,state in enumerate(NEB.states):
				if i==0 or i==len(NEB.states)-1: continue #don't change first or last state
				for j,b in enumerate(state):
					if j in spring_atoms:
						a,c = NEB.states[i-1][j], NEB.states[i+1][j]
						
						#find tangent
						tplus = np.array( [ c.x-b.x, c.y-b.y, c.z-b.z ] )
						tminus = np.array( [ a.x-b.x, a.y-b.y, a.z-b.z ] )
						dVmin = min(V[i+1]-V[i], V[i-1]-V[i])
						dVmax = max(V[i+1]-V[i], V[i-1]-V[i])
						if V[i+1]>V[i] and V[i]>V[i-1]: #not at an extremum, trend of V is up
							tangent = tplus
						elif V[i+1]<V[i] and V[i]<V[i-1]: #not at an extremum, trend of V is down
							tangent = tminus
						elif V[i+1]>V[i-1]: #at local extremum, next V is higher
							tangent = tplus*dVmax + tminus*dVmin
						else: #at local extremum, previous V is higher
							tangent = tplus*dVmin + tminus*dVmax
								
						#normalize tangent
						tangent /= np.linalg.norm(tangent)
						
						#find spring forces parallel to tangent
						F_spring_parallel = NEB.k*( utils.dist(c,b) - utils.dist(b,a) ) * tangent
						
						#find DFT forces perpendicular to tangent
						real_force = np.array( [b.fx,b.fz,b.fz] )
						F_real_perpendicular = real_force - np.dot(real_force, tangent)
						
						#set NEB forces
						b.fx, b.fy, b.fz = F_spring_parallel + F_real_perpendicular
						
						#TODO: get NEB energies rigorously such that NEB force = -gradient(NEB energy)
						#this is good enough for now:
						#energies[i] += 0.5*NEB.k*(utils.dist_squared(a,b) + utils.dist_squared(b,c))
						
						#old method
						#b.fx += NEB.k*(a.x-b.x) + NEB.k*(c.x-b.x)
						#b.fy += NEB.k*(a.y-b.y) + NEB.k*(c.y-b.y)
						#b.fz += NEB.k*(a.z-b.z) + NEB.k*(c.z-b.z)
						#energies[i] += 0.5*NEB.k*(utils.dist_squared(a,b) + utils.dist_squared(b,c))
			#set error
			NEB.error = sum(energies)
			#set forces
			NEB.forces = []
			for state in NEB.states[1:-1]:
				for a in state:
					NEB.forces += [-a.fx, -a.fy, -a.fz] #gradient of NEB.error
			#print data
			V = V[:1] + [ (e-V[0])/0.001 for e in V[1:] ]
			print NEB.step, '%8.6g:' % NEB.error, '%7.5g +' % V[0], ('%5.1f '*len(V[1:])) % tuple(V[1:])
			#increment step
			NEB.step += 1
	
		@staticmethod
		def get_error(coords):
			if NEB.error is None:
				NEB.calculate(coords)
			error = NEB.error
			NEB.error = None #set to None so it will recalculate next time
			return error
	
		@staticmethod
		def get_forces(coords):
			if NEB.forces is None:
				NEB.calculate(coords)
			forces = NEB.forces
			NEB.forces = None #set to None so it will recalculate next time
			return np.array(forces)*1.8 #convert from Hartree/Bohr to Hartree/Angstrom

	n = NEB(name, states, theory, k)
	minimize(NEB.get_error, np.array(NEB.coords_start), method='BFGS', jac=NEB.get_forces, options={'disp': True})


def neb_verlet(name, states, theory, extra_section='', queue=None, spring_atoms=None, k=0.1837, fit_rigid=True, dt=1.0): #Nudged Elastic Band. k for VASP is 5 eV/Angstrom, ie 0.1837 Hartree/Angstrom. 
#Cite NEB: http://scitation.aip.org/content/aip/journal/jcp/113/22/10.1063/1.1323224
	from scipy.optimize import minimize
	import numpy as np
	#set which atoms will be affected by virtual springs
	if not spring_atoms:#if not given, select all
		spring_atoms = range(len(states[0]))
	elif type(spring_atoms)==str: #a list of element names
		elements = spring_atoms.split()
		spring_atoms = [i for i,a in enumerate(states[0]) if a.element in elements]
	#class to contain working variables
			
	if fit_rigid:
		#center all states around spring-held atoms
		for s in states:
			center_x = sum([a.x for i,a in enumerate(s) if i in spring_atoms])/len(spring_atoms)
			center_y = sum([a.y for i,a in enumerate(s) if i in spring_atoms])/len(spring_atoms)
			center_z = sum([a.z for i,a in enumerate(s) if i in spring_atoms])/len(spring_atoms)
			for a in s:
				a.x -= center_x
				a.y -= center_y
				a.z -= center_z
		#rotate all states to be as similar to their neighbors as possible
		from scipy.linalg import orthogonal_procrustes
		for i in range(1,len(states)): #rotate all states to optimal alignment
			#only count spring-held atoms for finding alignment
			spring_atoms_1 = [(a.x,a.y,a.z) for j,a in enumerate(states[i]) if j in spring_atoms]
			spring_atoms_2 = [(a.x,a.y,a.z) for j,a in enumerate(states[i-1]) if j in spring_atoms]
			rotation = orthogonal_procrustes(spring_atoms_1,spring_atoms_2)[0]
			#rotate all atoms into alignment
			for a in states[i]:
				a.x,a.y,a.z = utils.matvec(rotation, (a.x,a.y,a.z))
	
	for s in states:
		for a in s:
			a.vx, a.vy, a.vz = 0.0, 0.0, 0.0
			a.ax, a.ay, a.az = 0.0, 0.0, 0.0
	for step in range(10):
		#start DFT jobs
		running_jobs = []
		for i,state in enumerate(states):
			if step>0:
				if (i==0 or i==len(states)-1): #only perform evaluations of endpoints on first step
					continue
				guess = ' Guess=Read'
			else: guess = '' #no previous guess for first step
			running_jobs.append( job('%s-%d-%d'%(name,step,i), theory+' Force'+guess, state, queue=queue, force=True, previous=('%s-%d-%d'%(name,step-1,i)) if step>0 else None, extra_section=extra_section) )
		#wait for jobs to finish
		for j in running_jobs: j.wait()
		#get forces and energies from DFT calculations
		V = [] #potential energy
		for i,state in enumerate(states):
			try:
				if (i==0 or i==len(states)-1): # if at an endpoint, just use first step's result
					step_to_use = 0
				else:
					step_to_use = step
				new_energy, new_atoms = parse_atoms('%s-%d-%d' % (name, step_to_use, i), check_convergence=False)
			except:
				print "Unexpected error in 'parse_atoms':", sys.exc_info()[0]
				print 'Job failed: %s-%d-%d'%(name,step,i); exit()
			V.append(new_energy)
			for a,b in zip(state, new_atoms):
				a.fx = b.fx; a.fy = b.fy; a.fz = b.fz
		#add spring forces to atoms
		for i,state in enumerate(states):
			if i==0 or i==len(states)-1: continue #don't change first or last state
			for j,b in enumerate(state):
				if j in spring_atoms:
					a,c = states[i-1][j], states[i+1][j]
					#find tangent
					tplus = np.array( [ c.x-b.x, c.y-b.y, c.z-b.z ] )
					tminus = np.array( [ a.x-b.x, a.y-b.y, a.z-b.z ] )
					dVmin = min(V[i+1]-V[i], V[i-1]-V[i])
					dVmax = max(V[i+1]-V[i], V[i-1]-V[i])
					if V[i+1]>V[i] and V[i]>V[i-1]: #not at an extremum, trend of V is up
						tangent = tplus
					elif V[i+1]<V[i] and V[i]<V[i-1]: #not at an extremum, trend of V is down
						tangent = tminus
					elif V[i+1]>V[i-1]: #at local extremum, next V is higher
						tangent = tplus*dVmax + tminus*dVmin
					else: #at local extremum, previous V is higher
						tangent = tplus*dVmin + tminus*dVmax
					#normalize tangent
					tangent /= np.linalg.norm(tangent)
					#find spring forces parallel to tangent
					F_spring_parallel = k*( utils.dist(c,b) - utils.dist(b,a) ) * tangent
					#find DFT forces perpendicular to tangent
					real_force = np.array( [b.fx,b.fz,b.fz] )
					F_real_perpendicular = real_force - np.dot(real_force, tangent)
					#set forces
					b.fx, b.fy, b.fz = F_spring_parallel + F_real_perpendicular
		# Projected Velocity Verlet, http://theory.cm.utexas.edu/henkelman/pubs/henkelman00_269.pdf
		masses = {'H':1.0008,'C':12.001,'O':15.998,'Si':28.085,'S':32.06,'Pb':208.2}
		for state in states:
			for a in state:
				#a = F/m
				a.ax_new = a.fx/masses[a.element]
				a.ay_new = a.fy/masses[a.element]
				a.az_new = a.fz/masses[a.element]
				#project velocity along direction of force
				a.vx, a.vy, a.vz = np.linalg.norm([a.vx,a.vy,a.vz]) * np.array([a.fx,a.fy,a.fz]) / np.linalg.norm([a.fx,a.fy,a.fz])
				#zero velocity if it points opposite to force
				if a.vx*a.fx < 0.0: a.vx = 0.0
				if a.vy*a.fy < 0.0: a.vy = 0.0
				if a.vz*a.fz < 0.0: a.vz = 0.0
				#update position
				a.x += a.vx*dt + 0.5*a.ax * dt**2
				a.y += a.vy*dt + 0.5*a.ay * dt**2
				a.z += a.vz*dt + 0.5*a.az * dt**2
				#update velocity
				a.vx += (a.ax + a.ax_new)*0.5*dt
				a.vy += (a.ay + a.ay_new)*0.5*dt
				a.vz += (a.az + a.az_new)*0.5*dt
				#update acceleration
				a.ax, a.ay, a.az = a.ax_new, a.ay_new, a.az_new
		#print data
		V = V[:1] + [ (e-V[0])/0.001 for e in V[1:] ]
		print step, '%7.5g +' % V[0], ('%5.1f '*len(V[1:])) % tuple(V[1:])
		#increment step
		step += 1


def optimize_pm6(name, examples, param_string, starting_params, queue=None): #optimize a custom PM6 semi-empirical method based on Gaussian examples at a higher level of theory
	from scipy.optimize import minimize
	import numpy as np
	
	examples = [utils.Struct(name=example, atoms=atoms(example)) for example in examples]
	for e in examples:
		e.bonds = [ (b.atoms[0].index-1, b.atoms[1].index-1) for b in utils.get_bonds(e.atoms) ]
	
	n_bonds = sum( [len(e.bonds) for e in examples] )
	
	counter = [0]
	def pm6_error(params):
		#run Gaussian jobs with new parameters
		for i,example in enumerate(examples):
			running_jobs = [ job('%s-%d-%d' % (name,counter[0],i), 'PM6=(Input,Print) Opt=Loose', example.atoms, extra_section=param_string%tuple(params), queue=queue, force=True)  ]
		#wait for all jobs to finish
		for j in running_jobs:
			j.wait()
		
		#get forces and energies resulting from new parameters
		geom_error = 0.0
		force_error = 0.0
		for i,example in enumerate(examples):
			try:
				new_energy, new_atoms = parse_atoms('%s-%d-%d'%(name,counter[0],i), check_convergence=False)
			except:
				print '%s-%d-%d'%(name,counter[0],i), 'has no data'
				exit()
			if parse_atoms('%s-%d-%d'%(name,counter[0],i)) is None:
				print '%s-%d-%d'%(name,counter[0],i), 'did not converge fully'
				#geom_error += 10.0 #discourage failure
		#compare results
			for b in example.bonds:
				d1 = utils.dist( example.atoms[b[0]], example.atoms[b[1]] )
				d2 = utils.dist( new_atoms[b[0]], new_atoms[b[1]] )
				geom_error += (d1-d2)**2
		
		error = geom_error/n_bonds
		
		print error**0.5, params
		
		counter[0]+=1
		
		return error
	
	minimize(pm6_error, starting_params, method='Nelder-Mead', options={'disp': True} )
	

