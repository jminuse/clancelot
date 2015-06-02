import os, re
import xml.etree.ElementTree as xml
import utils

def read_cml(name, parameter_file='oplsaa.prm', extra_parameters={}):
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
	for bond in root[1]:
		a, b = [int(n[1:]) for n in bond.attrib['atomRefs2'].split()]
		bonds.append( utils.Bond(atoms[a-1],atoms[b-1]) )
		atoms[a-1].bonded.append( atoms[b-1] )
		atoms[b-1].bonded.append( atoms[a-1] )
	
	angles, dihedrals = utils.get_angles_and_dihedrals(atoms)
	
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
	net_charge = sum([x.type.charge for x in atoms])
	if abs(net_charge)>0.01: print 'Non-neutral molecule, charge =', net_charge, ':', name; exit()
	
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
	

def write_cml(atoms, name=None):
	if name is None:
		name = 'out' #default filename is out.cml
	name += '.cml'
	
	#tree = xml.ElementTree()
	#tree._root = xml.fromstring('<molecule></molecule>')
	#tree.write(name)
	
	bonds = {}
	for a in atoms:
		for b in a.bonded:
			bonds[ tuple(sorted((a,b))) ] = True
	bonds = bonds.keys()
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
	atoms = []
	for line in open(name):
		columns = line.split()
		if len(columns)>=4:
			x,y,z = [float(s) for s in columns[1:4]]
			atoms.append( utils.Atom(element=columns[0], x=x, y=y, z=z, index=len(atoms)+1) )
	return atoms

def write_xyz(atoms, name_or_file=None):
	if not name_or_file:
		name_or_file = 'out' #default filename is out.xyz
	
	if type(name_or_file)==str: #if filename is provided, open it
		name = name_or_file
		f = open(name+'.xyz', 'w')
	else: #if open file is provided, just append to it
		f = name_or_file
	
	f.write(str(len(atoms))+'\nAtoms\n')
	for a in atoms:
		f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
	if name_or_file==name: #close file if we opened it
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


def write_lammps_data(system, pair_coeffs_included=True):
	atoms, bonds, angles, dihedrals, box_size, run_name = system.atoms, system.bonds, system.angles, system.dihedrals, system.box_size, system.name #unpack values from system

	dihedrals = [d for d in dihedrals if any(d.type.e)] #ignore dihedrals with no energy

	atom_types = dict( [(t.type,True) for t in atoms] ).keys() #unique set of atom types
	bond_types = dict( [(t.type,True) for t in bonds] ).keys()
	angle_types = dict( [(t.type,True) for t in angles] ).keys()
	dihedral_types = dict( [(t.type,True) for t in dihedrals] ).keys()
	
	atom_types.sort(key=lambda t:t.mass) # sort atom types by molecular weight
	
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] ) #assign unique numbers
	bond_type_numbers = dict( [(t,i+1) for i,t in enumerate(bond_types)] )
	angle_type_numbers = dict( [(t,i+1) for i,t in enumerate(angle_types)] )
	dihedral_type_numbers = dict( [(t,i+1) for i,t in enumerate(dihedral_types)] )
	
	f = open(run_name+'.data', 'w')
	f.write('LAMMPS Description\n\n%d atoms\n%d bonds\n%d angles\n%d dihedrals\n0  impropers\n\n' % (len(atoms), len(bonds), len(angles), len(dihedrals)) )
	f.write('%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n0  improper types\n' % (len(atom_types), len(bond_types), len(angle_types), len(dihedral_types)) )
	f.write('''
 -'''+str(box_size[0]/2)+''' '''+str(box_size[0]/2)+''' xlo xhi
 -'''+str(box_size[1]/2)+''' '''+str(box_size[1]/2)+''' ylo yhi
 -'''+str(box_size[2]/2)+''' '''+str(box_size[2]/2)+''' zlo zhi

Masses

'''+('\n'.join(["%d\t%f" % (atom_type_numbers[t], t.mass) for t in atom_types]))+'\n')

	if pair_coeffs_included: f.write('\nPair Coeffs\n\n'+('\n'.join(["%d\t%f\t%f" % (atom_type_numbers[t], t.vdw_e, t.vdw_r) for t in atom_types])) )

	if bonds: f.write("\n\nBond Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (bond_type_numbers[t], t.e, t.r) for t in bond_types]))
	if angles: f.write("\n\nAngle Coeffs\n\n"+'\n'.join(["%d\t%f\t%f" % (angle_type_numbers[t], t.e, t.angle) for t in angle_types]))
	if dihedrals: f.write("\n\nDihedral Coeffs\n\n"+'\n'.join(["%d\t%f\t%f\t%f\t%f" % ((dihedral_type_numbers[t],)+tuple(t.e)+((0.0,) if len(t.e)==3 else ()) ) for t in dihedral_types]))

	f.write("\n\nAtoms\n\n"+'\n'.join( ['\t'.join( [str(q) for q in [a.index, a.molecule_index, atom_type_numbers[a.type], a.type.charge, a.x, a.y, a.z]] ) for a in atoms]) ) #atom (molecule type charge x y z)

	if bonds: f.write('\n\nBonds\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, bond_type_numbers[b.type], b.atoms[0].index, b.atoms[1].index]]) for i,b in enumerate(bonds)]) ) #bond (type a b)
	if angles: f.write('\n\nAngles\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, angle_type_numbers[a.type]]+[atom.index for atom in a.atoms] ]) for i,a in enumerate(angles)]) ) #ID type atom1 atom2 atom3
	if dihedrals: f.write('\n\nDihedrals\n\n' + '\n'.join( ['\t'.join([str(q) for q in [i+1, dihedral_type_numbers[d.type]]+[atom.index for atom in d.atoms] ]) for i,d in enumerate(dihedrals)]) ) #ID type a b c d
	f.write('\n\n')
	f.close()


def packmol(system, molecules, molecule_ratio, density): #density in g/mL
	try:
		os.mkdir('packmol')
	except: pass
	os.chdir('packmol')
	
	f = open('pack.inp', 'w')
	f.write('''
	tolerance 2.0
	filetype xyz
	output out.xyz
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
	os.system('/fs/home/jms875/build/packmol/packmol < pack.inp')
	atoms = read_xyz('out.xyz')
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

