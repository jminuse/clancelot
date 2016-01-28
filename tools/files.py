import os, re
import xml.etree.ElementTree as xml
import utils

# Read the Chemical Markup Language (CML)
def read_cml(name, parameter_file='oplsaa.prm', extra_parameters={}, check_charges=True):
	if not name.endswith('.cml'):
		name += '.cml'
	tree = xml.parse(name)
	root = tree.getroot()
	
	atoms = []
	for atom in root[0]:
		a = utils.Atom( element=atom.attrib['elementType'], x=float(atom.attrib['x3']), y=float(atom.attrib['y3']), z=float(atom.attrib['z3']) )
		a.bonded = []
		a.index = int(atom.attrib['id'][1:])
		if 'formalCharge' in atom.attrib:
			a.type_index = int(atom.attrib['formalCharge'])
		if 'label' in atom.attrib:
			a.type_index = int(atom.attrib['label'])
		atoms.append(a)
	
	bonds = []
	if len(root)>1:
		for bond in root[1]:
			a, b = [int(n[1:]) for n in bond.attrib['atomRefs2'].split()]
			bonds.append( utils.Bond(atoms[a-1],atoms[b-1]) )
			atoms[a-1].bonded.append( atoms[b-1] )
			atoms[b-1].bonded.append( atoms[a-1] )
		angles, dihedrals = utils.get_angles_and_dihedrals(atoms)
	else:
		angles, dihedrals = [],[]
	
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
			for t in atom_types:
				if t.index==a.type_index:
					a.type = t
			for t in extra_parameters:
				if t==a.type_index:
					a.type = extra_parameters[t]
		
		#check charges
		if check_charges:
			net_charge = sum([x.type.charge for x in atoms])
			if abs(net_charge)>0.01: raise Exception('Non-neutral molecule, charge = %f: %s' % (net_charge, name))
		
		#set bond, angle, and dihedral types from parameter file
		for x in bonds+angles+dihedrals: 
			index2s = tuple([a.type.index2 for a in x.atoms])
			try:
				x.type = [t for t in bond_types+angle_types+dihedral_types if t.index2s==index2s or t.index2s==tuple(reversed(index2s))][0]
			except: pass

		#check consistency
		for x in bonds+angles+dihedrals: 
			if not x.type:
				print 'No type for structure indices', tuple([a.type.index2 for a in x.atoms]), ':', tuple([a.element for a in x.atoms]), ': atoms', tuple([a.index for a in x.atoms]),'in file', name
				if isinstance(x,utils.Dihedral): x.type = utils.Struct(index2s=index2s, e=(0.0,0.0,0.0)) #no params for dihedral is OK, just give warning
				else: print 'Exit'; exit()
	
	return atoms, bonds, angles, dihedrals
	

def write_cml(atoms, bonds, name=None):
	if name is None:
		name = 'out' #default filename is out.cml
	name += '.cml'
	
	#tree = xml.ElementTree()
	#tree._root = xml.fromstring('<molecule></molecule>')
	#tree.write(name)
	
	bonds = [ (b.atoms[0].index, b.atoms[1].index) for b in bonds ]
	bonds.sort()
	
	f = open(name, 'w')
	f.write('<molecule>\n <atomArray>\n')
	for i,a in enumerate(atoms):
		f.write('  <atom id="a%d" elementType="%s" x3="%f" y3="%f" z3="%f"' % (i+1, a.element, a.x, a.y, a.z) )
		try:
			f.write(' formalCharge="%d"' % a.type_index)
		except: pass
		try:
			f.write(' label="%d"' % a.label)
		except: pass
		f.write('/>\n')
	f.write(' </atomArray>\n <bondArray>\n')
	for b in bonds:
		f.write('  <bond atomRefs2="a%d a%d" order="1"/>\n' % b )
	f.write(' </bondArray>\n</molecule>')


def read_xyz(name):
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

def write_xyz(frames, name_or_file=None):
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
		f.write(str(len(atoms))+'\nAtoms\n')
		for a in atoms:
			f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
	if type(name_or_file)==str: #close file if we opened it
		f.close()

def read_opls_parameters(parameter_file):
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
	return elements, atom_types, bond_types, angle_types, dihedral_types


def write_lammps_data(system, pair_coeffs_included=False):
	#unpack values from system
	atoms, bonds, angles, dihedrals, box_size, run_name = system.atoms, system.bonds, system.angles, system.dihedrals, system.box_size, system.name
	#ignore dihedrals with no energy
	dihedrals = [d for d in dihedrals if any(d.type.e)]
	#get list of unique atom types
	atom_types = dict( [(t.type,True) for t in atoms] ).keys() #unique set of atom types
	bond_types = dict( [(t.type,True) for t in bonds] ).keys()
	angle_types = dict( [(t.type,True) for t in angles] ).keys()
	dihedral_types = dict( [(t.type,True) for t in dihedrals] ).keys()
	system.atom_types, system.bond_types, system.angle_types, system.dihedral_types = atom_types, bond_types, angle_types, dihedral_types
	# sort atom types by mass, largest masses first
	atom_types.sort(key=lambda t:-t.mass) 
	# get type numbers to identify types to LAMMPS
	for i,t in enumerate(atom_types): t.lammps_type = i+1
	for i,t in enumerate(bond_types): t.lammps_type = i+1
	for i,t in enumerate(angle_types): t.lammps_type = i+1
	for i,t in enumerate(dihedral_types): t.lammps_type = i+1
	#start writing file
	f = open(run_name+'.data', 'w')
	f.write('LAMMPS Description\n\n%d atoms\n%d bonds\n%d angles\n%d dihedrals\n0  impropers\n\n' % (len(atoms), len(bonds), len(angles), len(dihedrals)) )
	f.write('%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n0  improper types\n' % (len(atom_types), len(bond_types), len(angle_types), len(dihedral_types)) )
	f.write('''
 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

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


def packmol(system, molecules, molecule_ratio, density, seed=1): #density in g/mL
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
	''' % ((i,molecule_counts[i])+system.box_size) )
	f.close()
	os.system('/fs/home/jms875/build/packmol/packmol < '+system.name+'.packmol')
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

def inp_to_xyz(name, write=False,outName=None):
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
