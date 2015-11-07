import os, string, sys, re, shutil, copy
from subprocess import Popen
import utils, units, log, files

# Function that parses the route into the following situations
	# func/basis
	# key(a,b,c,...)
	# key(a)
	# key=a
	# key = a
	# key=(a,b,c,...)
def parse_route(route):
	parsed, parsed_list, r = [],[],route
	parsed += re.findall('(\w+\s*=\s*(?:(?:\([^\)]*\)|(?:\w+))))|(\w+\s*\([^\)]*\))|([A-z]\w+\/[A-z]\w+)',r)
	for p in parsed:
		for i in p:
			if i != '': parsed_list.append(i)
	for p in parsed_list: r = r.replace(p,'')
	return parsed_list, r.strip()

def chk_lint(run_name, route, atoms=[], extra_section='', queue='batch', procs=1, charge_and_multiplicity='0,1', title='run by gaussian.py', blurb=None, eRec=True, force=False, previous=None,neb=[False,None,None,None],err=False):
	# Check the extra section for required new lines and ****
	def chk_extra_section():
		# If no extra_section, don't do anything
		if extra_section == '': return True

		# Get the last section after splitting by '\n\n'.  Should be blank
		if extra_section.split('\n\n')[-1] == '': return True
		else:
			print('''Lint Err - Extra Section does not end in two new lines ('\\n\\n').''')
			return False

	def chk_pseudo_read(route, extra_section):
		# If we don't care about pseudo = read (required for Gen and GenECP), then ignore
		if route.split()[0].split('/')[1].lower() not in ['gen','genecp']: return True

		# Check if the route has pseudo=read
		r = route.strip().lower().replace(' ','')
		ss = ['pseudo=read','pseudo=(read']
		for s in ss:
			if r.find(s) > -1:
				has_pseudo=True
				break
			else: has_pseudo=False

		# Check if the extra section has the pseudo section (if it does, there should be 3: 1 for main mixed set, 2 for pseudo=read, 3 for end of section)
		if len(extra_section.split('\n\n'))==3: needs_pseudo=True
		else: needs_pseudo=False

		# Act accordingly
		if has_pseudo == needs_pseudo: return True
		else:
			print('Lint Err - Pseudo is incorrectly defined.')
			if has_pseudo: print('\tLinter found Psuedo=Read in route.')
			else: print('\tLinter did not find Pseudo=Read in route.')
			if needs_pseudo: print('\tLinter says Pseudo=Read required in route.')
			else: print('\tLinter says there should not be Pseudo=Read in route.')
			return False

	def chk_mixed_set(route, atoms, extra_section):
		# If we don't care about pseudo = read (required for Gen and GenECP), then ignore
		if route.split()[0].split('/')[1].lower() not in ['gen','genecp']: return True

		# Get list of elements
		#tmp = re.findall('((?:[A-Z][a-z]?\s*)+)\s0',extra_section)
		#elem_list = []
		#for t in tmp: elem_list += t.split()
		elem_list = [item for l in  ([ l.split() for l in re.findall('((?:[A-Z][a-e]?\s*)+)\s0',extra_section)]) for item in l]

		a_elem = []
		# Check if all atoms in elements
		for a in atoms:
			if hasattr(a,'element'):
				if len(units.elem_i2s(a.element))>2: e = units.elem_i2s(a.element)[:units.elem_i2s(a.element).find('-')]
				else: e = ''
				if (units.elem_i2s(a.element) not in elem_list) and (e not in elem_list):
					print('Lint Err - Element %s is not defined in mixed basis set.' % units.elem_i2s(a.element))
					return False
				if (e != '') and (e not in a_elem): a_elem.append(e)
				elif units.elem_i2s(a.element) not in a_elem: a_elem.append(units.elem_i2s(a.element))
			else:
				for b in a:
					if len(units.elem_i2s(b.element))>2: e = units.elem_i2s(b.element)[:units.elem_i2s(b.element).find('-')]
					else: e = ''
					if (units.elem_i2s(b.element) not in elem_list) and (e not in elem_list):
						print('Lint Err - Element %s is not defined in mixed basis set.' % units.elem_i2s(b.element))
						return False
					if (e != '') and (e not in a_elem): a_elem.append(e)
					elif units.elem_i2s(b.element) not in a_elem: a_elem.append(units.elem_i2s(b.element))	
		# Check if all atoms defined in mixed basis set are account for in elements
		for a in elem_list:
			if a not in a_elem:
				print('Lint Err - Element %s in mixed basis set is unaccounted for.' % a)
				return False

		return True
		
	funcs = [
	chk_extra_section, # Check if extra_section ends in '\n\n'
	chk_pseudo_read(route, extra_section), # Check if pseudo=Read error exists
	chk_mixed_set(route, atoms, extra_section) # Checks if all elements defined have appropriate basis sets
	# More functions placed here
	]

	if False in funcs:
		print('Errors were returned for %s.' % run_name)
		print('If you believe this/these error(s) to be wrong, please resubmit with lint=False.')
		sys.exit(0)

def job(run_name, route, atoms=[], extra_section='', queue='batch', procs=1, verbosity='N', charge_and_multiplicity='0,1', title='run by gaussian.py', blurb=None, eRec=True, force=False, previous=None, neb=[False,None,None,None], err=False, lint=False, mem=25):
	if lint: chk_lint(run_name, route, atoms, extra_section, queue, procs, charge_and_multiplicity, title, blurb, eRec, force, previous,neb,err)
	log.chk_gaussian(run_name,force=force,neb=neb) # Checks if run exists
	head = '#'+verbosity.upper()+' '+route+'\n\n'+title+'\n\n'+charge_and_multiplicity+'\n' # Header for inp file

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
%Mem=$$MEM$$MW
'''.replace("$$MEM$$",str(mem))
			inp.write(csh+head+xyz+extra_section+'\neof\nrm /tmp/*.rwf')
		if previous:	
			shutil.copyfile(previous+'.chk', run_name+'.chk')
		process_handle = Popen('/bin/csh %s.inp' % run_name, shell=True)
	pyPath = '../'+sys.argv[0]
	if not os.path.isfile(pyPath): pyPath = '../'+sys.argv[0][sys.argv[0].rfind('/')+1:] # This is if sys.argv[0] is a full path
	if not err: shutil.copyfile(pyPath, run_name+'.py')
	os.chdir('..')

	log.put_gaussian(run_name,route,extra_section=extra_section,blurb=blurb,eRec=eRec,force=force,neb=neb) # Places info in log if simulation is run
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

def parse_atoms(input_file, get_atoms=True, get_energy=True, check_convergence=True, get_time=False, counterpoise=False, parse_all=False):
	"""
	@input_file [str] : string name of log file

	Returns: (? energy, ? atoms, ? time) | None
	@energy [float] : If get_energy or parse_all, otherwise return omitted.
	@atoms |[atom list] : Iff parse_all, returns atom list list.
		   |[atom list list] : Iff not parse_all and get_atoms, atom list. Otherwise omitted.
	@time [float] : If get_time returns float (seconds). Otherwise, return omitted.

	Note that None may be returned in the event that Gaussian did not terminate normally (see 7 lines down).
	"""
	if input_file[-4:] != '.log':
		input_file = 'gaussian/'+input_file+'.log'
	contents = open(input_file).read()
	time = None	

	if check_convergence and get_energy and not parse_all and 'Normal termination of Gaussian 09' not in contents:
		return None
	if ('Normal termination of Gaussian 09' in contents) and (get_time | parse_all):
		m = re.search('Job cpu time: +(\S+) +days +(\S+) +hours +(\S+) +minutes +(\S+) +seconds', contents)
		try:
			time = float(m.group(1))*24*60*60 + float(m.group(2))*60*60 + float(m.group(3))*60 + float(m.group(4))
		except:
			pass

	if 'Summary of Optimized Potential Surface Scan' in contents and not parse_all:
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
		
	elif get_energy and not parse_all:
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

	if parse_all:
		energies = []
		atom_frames = []
		start=0
		orientation = 'Input orientation:'
		while True:
			try: #match energy
				input_orientation = contents.find(orientation, start)
				if input_orientation == -1:
					orientation = 'Standard orientation'
					print("\nWarning - No available Input Orientation, defaulting to Standard")
					input_orientation = contents.find(orientation, start)
				if input_orientation >= 0:
					start = input_orientation
				next_coordinates = contents.index('Coordinates (Angstroms)', start)
				start = contents.index('SCF Done', start)
				energies.append(float( re.search('SCF Done: +\S+ += +(\S+)', contents[start:]).group(1) ) )
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
		return energies, atom_frames, time

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
		last_coordinates = contents.rindex('Coordinates (Angstroms)')
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

def bandgap(input_file, parse_all=True):
	if input_file[-4:] != '.log':
		input_file = 'gaussian/'+input_file+'.log'
	hold = open(input_file).read()

	bandgap = []
	s = 'The electronic state is'
	while hold.find(s) != -1:
		hold = hold[hold.find(s) + len(s):]
		tmp = hold[:hold.find('Condensed')].split('\n')[1:-1]
		for i,t in enumerate(tmp):
			t = t.split()
			if t[1] == 'virt.':
				if i == 0:
					print("Error in calculating bandgap. Lowest eigenvalue energy is empty.")
					sys.exit()
				bandgap.append(float(t[4]) - float(tp[-1]))
				break
			tp = t
		hold = hold[hold.find('\n'):]
	if not parse_all: bandgap = bandgap[-1]
	
	return bandgap

def convergence(input_file):
	s = open('gaussian/'+input_file+'.log').read()
	s = s[s.rfind("Converged?"):].split('\n')[1:5]
	convergence = []
	if s == ['']: return None
	for c in s:
		c, tmp = c.split(), []
		tmp.append(' '.join(c[0:2]))
		tmp.append(float(c[2]))
		tmp.append(float(c[3]))
		tmp.append(c[4])
		convergence.append(tmp)
		
	return convergence

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
	if not input_file.startswith('gaussian/'):
		input_file = 'gaussian/' + input_file + '.log'
	with open(input_file) as inp:
		contents = inp.read()

	if 'Normal termination of Gaussian 09' not in contents:
		return None
	
	if contents.find('Fitting point charges to electrostatic potential') == -1:
		charges = None
	else:
		start = contents.rindex('Fitting point charges to electrostatic potential')
		end = contents.index('-----------------', start)
		charges = []
		for line in contents[start:end].splitlines():
			columns = line.split()
			if len(columns)==3:
				charges.append( float(columns[2]) )
	return charges

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
	
# A function that returns the binding energy of a molecule A with BSSE (Basis Set Superposition Error) corrections.
# job_total - This is the name of a gaussian job that holds the full system (optimized)
# job_A - This is the name of a gaussian job that holds the optimized molecule A
# job_B - This is the name of a gaussian job that holds the optimized molecule B
# zero_indexed_atom_indices_A - This is a list of indices for molecule A in job_total.  First values of a .xyz file start at 0.
def binding_energy(job_total, job_A, job_B, zero_indexed_atom_indices_A, route='SP SCRF(Solvent=Toluene)', blurb=None, procs=1, queue='batch', force=False, lint=False, bind_tck_name=None):
	#--------------------------
	def ghost_job(atoms, name, previous_job=None, route=None, blurb=None, procs=1, queue='batch', extras='', force=False, lint=False):
		# To ensure we do not overwrite a file we increment a value until we find that the run doesn't exist
		if os.path.isfile('gaussian/%s.inp' % name): # Check if normal file exists
			i = 1
			try: # Check if the last thing is an integer to increment
				int(name[name.rfind('_')+1:]) # If the last thing after _ is an integer, increment it
				while os.path.isfile('gaussian/%s_%d.inp' % (name[:name.rfind('_')],i)): i+=1
				name = name[:name.rfind('_')]+'_'+str(i)
			except ValueError:
				while os.path.isfile('gaussian/%s_%d.inp' % (name,i)): i+=1
				name = '%s_%d' % (name,i)		

		# Get appropriate level of theory, if not supplied
		if route is None:
			theory = open('gaussian/'+previous_job+'.inp').readline()[2:].strip().split()[0].split('/')
			# If we need to get mixed basis sets
			if theory[1].lower() in ['genecp','gen']:
				extras = open('gaussian/'+previous_job+'.inp').read().split('\n\n')[3:]
				extras = '\n\n'.join(extras)
				extras = extras.strip() + '\n\n'
				if len(extras.split('\n\n')) == 3: route = 'Pseudo=Read ' + route
			route = '/'.join(theory) + ' ' + route
		

		# Run the job and return the job name for the user to use later
		job(name, route, atoms=atoms, queue=queue, extra_section=extras, blurb=blurb, procs=procs, previous=previous_job, force=force, lint=lint)

		return name
	#--------------------------
	AB = atoms(job_total) # First get the atoms from the gaussian job for the full system
	AB_A = copy.deepcopy(AB)
	for i,atom in enumerate(AB_A): # For AB_A, we want all atoms not part of molecule A to be ghost atoms
		if i not in zero_indexed_atom_indices_A: atom.element+='-Bq'
	AB_B = copy.deepcopy(AB)
	for i,atom in enumerate(AB_B): # For AB_B we want all atoms part of molecule A to be ghost atoms
		if i in zero_indexed_atom_indices_A: atom.element+='-Bq'
	
	# Now AB_A is A from AB, AB_B is B from AB
	name1 = ghost_job(AB_A, job_total + '_A0', blurb=blurb, queue=queue, procs=procs, previous_job=job_total, force=force, lint=lint, route=route)
	name2 = ghost_job(AB_B, job_total + '_B0', blurb=blurb, queue=queue, procs=procs, previous_job=job_total, force=force, lint=lint, route=route)
	
	#non-rigid correction:
	AB_A = [atom for atom in AB_A if not atom.element.endswith('-Bq')]
	AB_B = [atom for atom in AB_B if not atom.element.endswith('-Bq')]
	name3 = ghost_job(AB_A, job_A + '_AB0', blurb=blurb, queue=queue, procs=procs, previous_job=job_A, force=force, lint=lint, route=route)
	name4 = ghost_job(AB_B, job_B + '_AB0'+('2' if job_A==job_B else '') , blurb=blurb, queue=queue, procs=procs, previous_job=job_B, force=force, lint=lint, route=route)
	
	# To get the binding energy we need to take into account the superposition error and the deformation error:
	# Superposition Error Correction is done by taking the total energy of the job and subtracting from it:
	#   	name1 = Original system with molecule A as ghost atoms and basis set of original system
	#		name2 = Original system with molecule B as ghost atoms and basis set of original system
	# Deformation energy corrections are done by taking the difference in energy of molecules in the same basis set between
	# geometries:
	#		Add the following
	#		name3 = Molecule A in the geometry of the original System but in the basis set of what molecule A was done in
	#		name4 = Molecule B in the geometry of the original System but in the basis set of what molecule B was done in
	# 		Subtract the following
	#		job_A = Molecule A in the basis set it was done in
	#		job_B = Molecule B in the basis set it was done in
	print 'E_binding = %s - %s - %s + %s + %s - %s - %s' % (job_total, name1, name2, name3, name4, job_A, job_B)

	if bind_tck_name==None:
		i = 0
		while os.path.isfile('bind_tck_#.py'.replace('#',str(i))): i += 1
		
		f = open('bind_tck_#.py'.replace('#',str(i)),'w')
	else:
		if len(bind_tck_name) > 3 and bind_tck_name[-3:] == '.py':
			f = open(bind_tck_name,'w')
		else:
			f = open("%s.py" % bind_tck_name,'w')

	job_names = [job_total, name1, name2, name3, name4, job_A, job_B]
	for i in range(len(job_names)): job_names[i] = "'"+job_names[i]+"'"
	job_names = ', '.join(job_names)

	f.write('''from merlin import *
try:
	s_units=sys.argv[1]
except:
	s_units='Ha'
job_names = [$$$$$]
# Check if jobs are still running
for s in log.get_jlist():
	if s in job_names:
		print("Sorry, all simulations haven't finished yet...")
		sys.exit()

# Else, we can get the energy
energies = []
print('Jobs Calculated From: ')
for s in job_names:
	try:
		e,atoms = g09.parse_atoms(s)
		energies.append(e)
		print '\t%20s\t%lg' % (s,e), ' '.join([a.element for a in atoms])
	except:
		raise Exception('Error, could not get data from %s.' % s)

sp_corr = units.convert_energy('Ha',s_units,energies[0] - energies[1] - energies[2])
deform_a = units.convert_energy('Ha',s_units,energies[3] - energies[5])
deform_b = units.convert_energy('Ha',s_units,energies[4] - energies[6])
geom_corr = deform_a + deform_b
print('------------')
print('Superposition Correction = '+str(sp_corr)+' '+s_units)
print('Geometry Correction = '+str(geom_corr)+' '+s_units)
print('\tDeformation of A = '+str(deform_a)+' '+s_units)
print('\tDeformation of B = '+str(deform_b)+' '+s_units)
print('Binding Energy = '+str(sp_corr + geom_corr)+' '+s_units)'''.replace('$$$$$',job_names))

	f.close()

def read(input_file):
	data = utils.DFT_out(input_file, 'g09')

	data.frames = parse_atoms(input_file, get_atoms=True, get_energy=False, check_convergence=False, get_time=False, counterpoise=False, parse_all=True)[1]
	if data.frames == []: data.frames = None
	data.atoms = data.frames[-1] if type(data.frames)==list and type(data.frames[0])==list else data.frames
	data.energies = parse_atoms(input_file, get_atoms=False, get_energy=True, check_convergence=False, get_time=False, counterpoise=False, parse_all=True)[0]
	data.charges_CHELPG = parse_chelpg(input_file)
	data.charges = data.charges_CHELPG
	data.convergence = convergence(input_file)
	data.converged = parse_atoms(input_file, get_atoms=False, get_energy=True, check_convergence=True, get_time=False, counterpoise=False, parse_all=False) is not None
	data.time = parse_atoms(input_file, get_atoms=False, get_energy=False, check_convergence=False, get_time=True, counterpoise=False, parse_all=True)[2]
	data.bandgap = bandgap(input_file)

	return data