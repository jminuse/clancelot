import os, re
import xml.etree.ElementTree as xml
import utils
import datetime
from warnings import warn

# Read the Chemical Markup Language (CML)
def read_cml(name, parameter_file='oplsaa.prm', extra_parameters={}, test_charges=True, allow_errors=False):
	if not name.endswith('.cml'):
		name += '.cml'
	tree = xml.parse(name)
	root = tree.getroot()

	# If periodic box information was passed to cml file, first root is the periodic box
	# If no periodic box info, first root is atomArray
	# Fill appropriate section when it is found
	atoms, bonds, angles, dihedrals = [], [], [], []

	for child in root:
		# Skip crystal information
		if child.tag == 'crystal':
			continue

		# Add atoms
		if child.tag == 'atomArray':
			for atom in child:
				a = utils.Atom( element=atom.attrib['elementType'], x=float(atom.attrib['x3']), y=float(atom.attrib['y3']), z=float(atom.attrib['z3']) )
				a.bonded = []
				a.index = int(atom.attrib['id'][1:])
				#if 'formalCharge' in atom.attrib:
				#	a.type_index = int(atom.attrib['formalCharge'])
				if 'label' in atom.attrib:
					a.type_index = int(atom.attrib['label'])
				atoms.append(a)

		# Add atoms
		if child.tag == 'bondArray':
			for bond in child:
				a, b = [int(n[1:]) for n in bond.attrib['atomRefs2'].split()]
				bonds.append( utils.Bond(atoms[a-1],atoms[b-1]) )
				atoms[a-1].bonded.append( atoms[b-1] )
				atoms[b-1].bonded.append( atoms[a-1] )
			angles, dihedrals = utils.get_angles_and_dihedrals(atoms)

	if parameter_file:
		atoms, bonds, angles, dihedrals = set_forcefield_parameters(atoms, bonds=bonds, angles=angles, dihedrals=dihedrals, name=name, parameter_file=parameter_file, extra_parameters=extra_parameters, test_charges=test_charges, allow_errors=allow_errors)

	return atoms, bonds, angles, dihedrals

# 1st param: either a list of util.Atom objects, a utils.Molecule object or a utils.System object
def write_cml(atoms_or_molecule_or_system, bonds=[], name=None):
	if name is None:
		name = 'out' #default filename is out.cml
	if not name.endswith('.cml'):
		name += '.cml'

	# Determine whether input is an atom list or a system object
	# If it is a system object, compile information to write cml file
	if isinstance(atoms_or_molecule_or_system, utils.System):
		system = atoms_or_molecule_or_system
		atoms, bonds, periodic = system.atoms, system.bonds, system.periodic
	elif isinstance(atoms_or_molecule_or_system, utils.Molecule):
		molecule = atoms_or_molecule_or_system
		atoms = molecule.atoms
		bonds = bonds or molecule.bonds
		periodic = False
	elif isinstance(atoms_or_molecule_or_system[0], utils.Atom):
		atoms =  atoms_or_molecule_or_system
		periodic = False
	else:
		raise Exception('Unable to write cml file = %s' % (name))

	for i,a in enumerate(atoms):
		a.index_from_1 = i+1
	bond_indices = [ (b.atoms[0].index_from_1 , b.atoms[1].index_from_1) for b in bonds ]
	bond_indices.sort()

	f = open(name, 'w')
	f.write('<molecule>\n')

	# Write crystal information if provided. Assign it is a crystal if it is periodic
	if periodic:
		f.write(' <crystal>\n')
		f.write('  <scalar title="a" units="units:angstrom">%3.6f</scalar>\n' % (system.box_size[0]) )
		f.write('  <scalar title="b" units="units:angstrom">%3.6f</scalar>\n' % (system.box_size[1]) )
		f.write('  <scalar title="c" units="units:angstrom">%3.6f</scalar>\n' % (system.box_size[2]) )
		f.write('  <scalar title="alpha" units="units:degree">%3.6f</scalar>\n' % (system.box_angles[0]) )
		f.write('  <scalar title="beta" units="units:degree">%3.6f</scalar>\n' % (system.box_angles[1]) )
		f.write('  <scalar title="gamma" units="units:degree">%3.6f</scalar>\n' % (system.box_angles[2]) )
		f.write(' </crystal>\n')

	# Write atom information
	f.write(' <atomArray>\n')
	for i,a in enumerate(atoms):
		f.write('  <atom id="a%d" elementType="%s" x3="%f" y3="%f" z3="%f"' % (i+1, a.element, a.x, a.y, a.z) )
		if hasattr(a, 'type.charge'):
			f.write(' formalCharge="%d"' % a.type.charge)

		if hasattr(a, 'type_index') and a.type_index!=None:
			f.write(' label="%d"' % a.type_index)
		elif hasattr(a, 'label'):
			f.write(' label="%s"' % str(a.label))
		f.write('/>\n')
	f.write(' </atomArray>\n <bondArray>\n')
	for pair in bond_indices:
		f.write('  <bond atomRefs2="a%d a%d" order="1"/>\n' % pair )
	f.write(' </bondArray>\n</molecule>')

# Imports the atom style dump file from lammps. VMD automatically reads this when labelled as .lammpstrj
# Default is to import everything but you will get better performance if you turn off the data you do not need
def read_lammpstrj(name, read_atoms=True, read_timesteps=True, read_num_atoms=True, read_box_bounds=True):
	if not name.endswith('.lammpstrj') and '.' not in name:
		name += '.lammpstrj'

	# If file does not exist, return empty lammpstrj object
	if not os.path.isfile(name):
		warn('Expected lammps trajectory file does not exist at %s' % (name))
		data = ''
	else:
		data = open(name,'r').read()

	# Get all the positions
	section, frames = data, []
	s = 'ITEM: ATOMS id type x y z'
	while read_atoms and (s in section):
		section = section[section.find(s)+len(s):]
		atom_block = section[:section.find('\nITEM: TIMESTEP')].split('\n')[1:]
		frame = []
		for line in atom_block:
			a = line.split()

			# Check if atom has expected number of characteristics
			if len(a) == 5:
				frame.append(utils.Atom(a[1],float(a[2]),float(a[3]),float(a[4]),index=a[0]))
			else:
				print('Atom skipped due to missing information')

		frames.append(frame)

	if frames:
		atoms = frames[-1]
	else:
		atoms = None

	# Get all timesteps
	section, timesteps = data, []
	s = 'ITEM: TIMESTEP'
	while read_timesteps and (s in section):
		section = section[section.find(s)+len(s):]
		tmp = section[:section.find('\nITEM: NUMBER OF ATOMS')].split('\n')[1:]
		for line in tmp:
			a = line.split()
			timesteps.append(int(a[0]))

	if len(timesteps) > 0:
		final_timestep = timesteps[-1]
	else:
		final_timestep = None

	# Get number of atoms. Useful if number of atoms change during simulation, such as during a deposition
	section, atom_counts = data, []
	s = 'ITEM: NUMBER OF ATOMS'
	while read_num_atoms and (s in section):
		section = section[section.find(s)+len(s):]
		tmp = section[:section.find('\nITEM: BOX BOUNDS')].split('\n')[1:]
		for line in tmp:
			a = line.split()
			atom_counts.append(int(a[0]))

	if len(atom_counts) > 0:
		atom_count = atom_counts[-1]
	else:
		atom_count = None

	# Get box bounds
	# Currently only imports orthorhombic crystal information aka all angles = 90 degrees
	section, box_bounds_list = data, []
	s = 'ITEM: BOX BOUNDS'
	while read_box_bounds and (s in section):
		section = section[section.find(s)+len(s):]
		tmp = section[:section.find('\nITEM: ATOMS')].split('\n')[1:]
		box_bounds = utils.Struct(xlo=None, xhi=None, ylo=None, yhi=None, zlo=None, zhi=None)

		for line in tmp:
			a = line.split()
			if box_bounds.xlo is None:
				box_bounds.xlo = float(a[0])
				box_bounds.xhi = float(a[1])
			elif box_bounds.ylo is None:
				box_bounds.ylo = float(a[0])
				box_bounds.yhi = float(a[1])
			elif box_bounds.zlo is None:
				box_bounds.zlo = float(a[0])
				box_bounds.zhi = float(a[1])

		box_bounds_list.append(box_bounds)

	if len(box_bounds_list) > 0:
		box_bounds = box_bounds_list[-1]
	else:
		box_bounds = None

	# Create object to store all results
	data = utils.sim_out(name, 'lammps')

	# Record all lammps trajectory data into results object
	data.frames = frames
	data.atoms = atoms
	data.timesteps = timesteps
	data.final_timestep = final_timestep
	data.atom_counts = atom_counts
	data.atom_count = atom_count
	data.box_bounds_list = box_bounds_list
	data.box_bounds = box_bounds
	data.last_modified = 'Null' # Stores when lammpstrj was last modified in seconds

	return data

def read_xyz(name):
	if not name.endswith('.xyz') and '.' not in name:
		name += '.xyz'
	lines = open(name).readlines()
	atom_count = int(lines[0].split()[0])
	lines_by_frame = [ lines[i:i+atom_count+2] for i in range(0,len(lines),atom_count+2) ]

	frames = []
	for frame in lines_by_frame:
		atoms = []
		for line in frame[2:]:
			columns = line.split()
			if len(columns)>=4:
				x,y,z = [float(s) for s in columns[1:4]]
				atoms.append( utils.Atom(element=columns[0], x=x, y=y, z=z, index=len(atoms)+1) )
		if len(atoms)>0: frames.append(atoms)

	if len(frames)==1:
		return frames[0]
	else:
		return frames

def write_xyz(frames_or_system, name_or_file=None, ID='Atoms'):
	# Determine whether input is an atom list (frames) or a system object
	# If it is a system object, compile information to write xyz file
	if isinstance(frames_or_system, utils.System):
		system = frames_or_system
		frames = system.atoms
	else:
		frames =  frames_or_system

	if not name_or_file:
		name_or_file = 'out' #default filename is out.xyz

	if type(name_or_file)==str: #if filename is provided, open it
		name = name_or_file
		f = open(name+'.xyz', 'w')
	else: #if open file is provided, just append to it
		f = name_or_file

	if isinstance(frames[0],utils.Atom):
		frames = [frames] #we want to write a list of frames, so make it one

	for atoms in frames:
		f.write(str(len(atoms))+'\n'+ID+'\n')
		for a in atoms:
			f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )

	if type(name_or_file)==str: #close file if we opened it
		f.close()

def read_opls_parameters(parameter_file='oplsaa.prm'):
	if read_opls_parameters.atom_types is None:
		elements = {}; atom_types = []; bond_types = []; angle_types = []; dihedral_types = []
		for line in open(parameter_file):
			columns = line.split()
			if not columns: continue
			if columns[0]=='atom':
				m = re.match('atom +(\d+) +(\d+) +(\S+) +"([^"]+)" +(\d+) +(\S+) +(\d+)', line)
				atom_type = utils.Struct(index=int(m.group(1)), index2=int(m.group(2)), element_name=m.group(3), notes=m.group(4), element=int(m.group(5)), mass=float(m.group(6)), bond_count=int(m.group(7) ) )
				if atom_type.element not in elements: elements[atom_type.element] = []
				elements[atom_type.element].append(atom_type)
				if '(UA)' in atom_type.notes:
					atom_type.element = 0 #reject united-atom parameters
				atom_types.append(atom_type)
			elif columns[0]=='vdw':
				atom_types[int(columns[1])-1].vdw_r = max( float(columns[2]), 1.0)
				atom_types[int(columns[1])-1].vdw_e = max( float(columns[3]), 0.01)
			elif columns[0]=='charge':
				atom_types[int(columns[1])-1].charge = float(columns[2])
			elif columns[0]=='bond':
				bond_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:3]]),e=float(columns[3]),r=float(columns[4])) )
			elif columns[0]=='angle':
				angle_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:4]]),e=float(columns[4]),angle=float(columns[5])) )
			elif columns[0]=='torsion':
				dihedral_types.append( utils.Struct(index2s=tuple([int(s) for s in columns[1:5]]),e=tuple([float(s) for s in columns[5::3]])) )
				if len(dihedral_types[-1].e)==3:
					dihedral_types[-1].e = dihedral_types[-1].e + (0.,)
		read_opls_parameters.elements = elements
		read_opls_parameters.atom_types = atom_types
		read_opls_parameters.bond_types = bond_types
		read_opls_parameters.angle_types = angle_types
		read_opls_parameters.dihedral_types = dihedral_types

	return read_opls_parameters.elements, read_opls_parameters.atom_types, read_opls_parameters.bond_types, read_opls_parameters.angle_types, read_opls_parameters.dihedral_types
#store these values between calls:
read_opls_parameters.elements = None
read_opls_parameters.atom_types = None
read_opls_parameters.bond_types = None
read_opls_parameters.angle_types = None
read_opls_parameters.dihedral_types = None

def set_forcefield_parameters(atoms, bonds=[], angles=[], dihedrals=[], parameter_file='oplsaa.prm', name='unnamed', extra_parameters={}, test_consistency=True, test_charges=True, allow_errors=False):
	elements, atom_types, bond_types, angle_types, dihedral_types = [], [], [], [], []

	# If parameter_file=None, only extra parameters will be passed.
	if parameter_file:
		elements, atom_types, bond_types, angle_types, dihedral_types = read_opls_parameters(parameter_file)

	#add extra parameters, if any
	for index2s,params in extra_parameters.items():
		if type(index2s)==int: continue #skip these
		if len(index2s)==2:
			bond_types.append( utils.Struct(index2s=index2s, e=params[0], r=params[1]) )
		elif len(index2s)==3:
			angle_types.append( utils.Struct(index2s=index2s, e=params[0], angle=params[1]) )
		elif len(index2s)==4:
			dihedral_types.append( utils.Struct(index2s=index2s, e=tuple(params)) )

	#set atom types
	for a in atoms:
		if a.type_index is None:
			raise Exception('OPLS label is missing from atom %d' % (a.index))
		#print(a.printSelf())
		for t in atom_types:
			if t.index==a.type_index:
				a.type = t
		for t in extra_parameters:
			if t==a.type_index:
				a.type = extra_parameters[t]

	#set bond, angle, and dihedral types from parameter file
	for x in bonds+angles+dihedrals:
		index2s = tuple([a.type.index2 for a in x.atoms])
		try:
			# Match type from opls parameters. Updated to allow for wildcard torsions
			for t in (bond_types+angle_types+dihedral_types):
				# Check for wildcard torsions
				if len(index2s) == 4 and len(t.index2s) == 4:
					if t.index2s[0] == 0 and t.index2s[3] == 0:
						match = t.index2s[1:3]==index2s[1:3] or t.index2s[1:3]==tuple(reversed(index2s))[1:3]
					elif t.index2s[0] == 0:
						match = t.index2s[1:4]==index2s[1:4] or t.index2s[1:4]==tuple(reversed(index2s))[1:4]
					elif t.index2s[3] == 0:
						match = t.index2s[0:3]==index2s[0:3] or t.index2s[0:3]==tuple(reversed(index2s))[0:3]
					else:
						match = t.index2s==index2s or t.index2s==tuple(reversed(index2s))

				# Check bonds and angles
				else:
					match = t.index2s==index2s or t.index2s==tuple(reversed(index2s))

				if match:
					x.type = t
					break

			#print([t for t in bond_types+angle_types+dihedral_types if t.index2s==index2s or t.index2s==tuple(reversed(index2s))])
			#print([t for t in bond_types+angle_types+dihedral_types if t.index2s==index2s or t.index2s==tuple(reversed(index2s))][0])
			
			#x.type = [t for t in bond_types+angle_types+dihedral_types if t.index2s==index2s or t.index2s==tuple(reversed(index2s))][0]
		except: pass

	if test_charges:
		check_net_charge(atoms, name=name)

	if test_consistency:
		check_consistency(atoms, bonds, angles, dihedrals, name=name, allow_errors=allow_errors)

	return atoms, bonds, angles, dihedrals

# Check charges to see if it is a neutral molecule
# Requires force field parameters: type.charge
# Raises exception if non-neutral. Consider changing to warning to allow non-neutral molecules
def check_net_charge(atoms, name=''):
	net_charge = sum([x.type.charge for x in atoms])
	if abs(net_charge)>0.01: raise Exception('Non-neutral molecule, charge = %f: %s' % (net_charge, name))

# Check to see if all possible force field parameters have been assigned.
# Raises exception if missing an bond or angle. Missing dihedrals allowed.
# Can turn off raising exceptions.
def check_consistency(atoms, bonds, angles, dihedrals, name='', allow_errors=False):
	for x in bonds+angles+dihedrals:
		# Compile all index types?
		index2s = tuple([a.type.index2 for a in x.atoms])

		if not x.type:
			print 'No type for structure indices', tuple([a.type.index2 for a in x.atoms]), ':', tuple([a.element for a in x.atoms]), ': atoms', tuple([a.index for a in x.atoms]),'in file', name
			if isinstance(x,utils.Dihedral): x.type = utils.Struct(index2s=index2s, e=(0.0,0.0,0.0)) #no params for dihedral is OK, just give warning
			elif allow_errors:
				continue
			else: print 'Exit'; exit()

#@profile
def write_lammps_data(system, pair_coeffs_included=False):
	#unpack values from system
	atoms, bonds, angles, dihedrals, box_size, box_angles, run_name = system.atoms, system.bonds, system.angles, system.dihedrals, system.box_size, system.box_angles, system.name
	#unpack lammps box parameters from system
	xlo, ylo, zlo, xhi, yhi, zhi, xy, xz, yz = system.xlo, system.ylo, system.zlo, system.xhi, system.yhi, system.zhi, system.xy, system.xz, system.yz

	#ignore angles with no energy
	angles = [ang for ang in angles if ang.type.e > 0]
	#ignore dihedrals with no energy
	dihedrals = [d for d in dihedrals if any(d.type.e)]

	#get list of unique atom types
	atom_types = dict( [(t.type,True) for t in atoms] ).keys() #unique set of atom types
	bond_types = dict( [(t.type,True) for t in bonds] ).keys()
	angle_types = dict( [(t.type,True) for t in angles] ).keys()
	dihedral_types = dict( [(t.type,True) for t in dihedrals] ).keys()
	system.atom_types, system.bond_types, system.angle_types, system.dihedral_types = atom_types, bond_types, angle_types, dihedral_types
	# sort atom types by mass, largest masses first
	atom_types.sort(key=lambda t:-t.mass + (-1e6 if hasattr(t,'reax') else 0) )
	# get type numbers to identify types to LAMMPS
	for i,t in enumerate(atom_types): t.lammps_type = i+1
	for i,t in enumerate(bond_types): t.lammps_type = i+1
	for i,t in enumerate(angle_types):
		#print(angle_types)
		#print(t)
		#print(t.lammps_type)
		t.lammps_type = i+1
	for i,t in enumerate(dihedral_types): t.lammps_type = i+1
	#start writing file
	f = open(run_name+'.data', 'w')
	f.write('LAMMPS Description\n\n%d atoms\n%d bonds\n%d angles\n%d dihedrals\n0  impropers\n\n' % (len(atoms), len(bonds), len(angles), len(dihedrals)) )
	f.write('%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n0  improper types\n' % (len(atom_types), len(bond_types), len(angle_types), len(dihedral_types)) )
	f.write('%3.5f %3.5f xlo xhi\n' % (xlo, xhi))
	f.write('%3.5f %3.5f ylo yhi\n' % (ylo, yhi))
	f.write('%3.5f %3.5f zlo zhi\n' % (zlo, zhi))
	# If the system is triclinic box
	if abs(box_angles[0]-90) > 0.001 or abs(box_angles[1]-90) > 0.001 or abs(box_angles[2]-90) > 0.001:
		f.write('%3.5f %3.5f %3.5f xy xz yz\n' % (xy, xz, yz))

	f.write('''
Masses

'''+('\n'.join(["%d\t%f" % (t.lammps_type, t.mass) for t in atom_types]))+'\n')

	if pair_coeffs_included: f.write('\nPair Coeffs\n\n'+('\n'.join(["%d\t%f\t%f" % (t.lammps_type, t.vdw_e, t.vdw_r) for t in atom_types])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (t.lammps_type, t.e, t.r) for t in bond_types]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (t.lammps_type, t.e, t.angle) for t in angle_types]))
	if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d\t%f\t%f\t%f\t%f" % ((t.lammps_type,)+tuple(t.e)+((0.0,) if len(t.e)==3 else ()) ) for t in dihedral_types]))

	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in [a.index, a.molecule_index, a.type.lammps_type, a.type.charge, a.x, a.y, a.z]] ) for a in atoms]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, b.type.lammps_type, b.atoms[0].index, b.atoms[1].index]]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, a.type.lammps_type]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, d.type.lammps_type]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()


def packmol(system, molecules, molecule_ratio=(1,), density=1.0, seed=1): #density in g/mL
	try:
		os.mkdir('packmol')
	except: pass
	os.chdir('packmol')

	f = open(system.name+'.packmol', 'w')
	f.write('''
	tolerance 2.0
	filetype xyz
	output '''+system.name+'''.packed.xyz
	seed '''+str(seed)+'''
	''')
	density *= 0.6022 # convert density to amu/angstrom^3. 1 g/mL = 0.6022 amu/angstrom^3
	average_molecular_weight = sum([a.type.mass*molecule_ratio[i] for i in range(len(molecules)) for a in molecules[i].atoms]) / sum(molecule_ratio)
	count = density * system.box_size[0] * system.box_size[1] * system.box_size[2] / average_molecular_weight
	molecule_counts = [int(round( count*x/sum(molecule_ratio) )) for x in molecule_ratio]
	for i,m in enumerate(molecules):

		xyz_file = open('%d.xyz' % i, 'w')
		xyz_file.write(str(len(m.atoms))+'\nAtoms\n')
		for a in m.atoms:
			xyz_file.write('%s%d %f %f %f\n' % (a.element, i, a.x, a.y, a.z) )
		xyz_file.close()

		f.write('''
	structure %d.xyz
	  number %d
	  inside box 0. 0. 0. %f %f %f
	end structure
	''' % ((i,molecule_counts[i])+tuple(system.box_size)) )
	f.close()
	os.system('/fs/home/jms875/build/packmol/packmol < '+system.name+'.packmol > packmol.log')
	atoms = read_xyz(system.name+'.packed.xyz')
	os.chdir('..')

	#now have a list of atoms with element = H0 for molecule 0, H1 for molecule 1, etc
	i = 0
	while i < len(atoms):
		molecule_number = int(atoms[i].element[-1])
		molecule = molecules[molecule_number]
		system.add(molecule)
		for a in system.atoms[-len(molecule.atoms):]: #update positions of latest molecule
			a.x, a.y, a.z = atoms[i].x, atoms[i].y, atoms[i].z
			i+=1

def inp_to_xyz(name, write=False, outName=None):
	warn(
        "this function is not used and will be removed soon.",
        DeprecationWarning
    )
	data = open("gaussian/"+name+".inp",'r').read().split('\n')

	# Get start of data
	for i,s in enumerate(data):
		try:
			if(s.split()[0]=='run'): break
		except:
			pass
	i += 3

	# Get end of data
	for j,s in enumerate(data[i:]):
		if s == '': break

	data = data[i:j+i]

	if write:
		if outName == None: f = open('inp_'+name+'.xyz','w')
		else: f = open(outName,'w')
		f.write(str(len(data))+'\n')
		f.write('Atoms'+'\n')
		for s in data: f.write(s+'\n')
		f.write('\n')
		f.close()

	atoms = []
	for i,d in enumerate(data):
		d = d.split()
		atoms.append(utils.Atom(element=d[0], x=float(d[1]), y=float(d[2]), z=float(d[3]), index=i))

	return atoms

# Returns the last time a file was modified in seconds
def last_modified(name):
	if not os.path.isfile(name):
		warn('Expected file does not exist at %s' % (name))
		return 0

	statinfo = os.stat(name)
	return datetime.datetime.fromtimestamp(statinfo.st_mtime)