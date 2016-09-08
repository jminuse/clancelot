import os, sys
import math, copy, subprocess, time, re
import numpy
import files, constants
import numpy as np
import scipy
from units import elem_i2s
from warnings import warn, simplefilter
from sysconst import opls_path
import random


simplefilter('always', DeprecationWarning)

class Struct(object):
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
	def __repr__(self):
		return str( dict([ (a,None) if type(self.__dict__[a]) in (list,dict) else (a,self.__dict__[a]) for a in self.__dict__]) )


class _Physical(object):
	"""
	An instance of Physical is a physical object or collection of physical objects,
	such as a Molecule, Atom, Bond, Dihedral, System, etc.
	This class is designed as a parent class to Atom, Molecule, Bond, System, and
	Dihedral to provide a basic equals() method and other shared methods.

	Invariant: Everything about a _Physical instance must be described by attributes
	whose order does not vary or whose attributes' equivalence is not defined by
	their order (e.g. Sets). No _Physical instance may have a 2-D list as an
	attribute, or equals() will fail.
	"""

	def equals(self, other):
		"""
		Basic equivalence of self and other if they have the same pointer id, or if all
		non-callable dir entries of the two objects are equivalent.

		Since several parameters are not preserved between intuitively equivalent
		_Physical instances (e.g. Atom.index or Atom.molecule_index, Angle.theta
		because of how Angles are added to a System), and these must be explicitly
		discounted in this method.

		The current list of explicitly discounted parameters are:
		Angle.theta
		Atom.index
		Atom.molecule_index
		Atom.bonded's Atom identities: the number of atoms in the list is
		compared instead.
		"""
		#Check basic identifiers to see if the two are obviously equal or inequal
		if id(self) == id(other):
			return True
		if type(self) != type(other):
			return False

		#Compare all attributes of the two objects. If they are exactly identical,
		#Return true.
		otherDict = other.__dict__
		selfDict = self.__dict__

		if selfDict == otherDict:
			return True

		#If the two objects don't have the same attribute names, return false
		if set(selfDict.keys()) != set(otherDict.keys()):
			return False

		#Go through each attribute individually, and check equality. If they
		#are _Physical instances, use .equals() to compare.
		for a in selfDict:

			#Simple attribute checking
			if selfDict[a] == otherDict[a]:
				continue

			#Passed-on attributes which are not conserved between atoms
			elif a == 'theta' and isinstance(self,Angle):
				continue
			elif isinstance(self,Atom) and (a=='index' or a=='molecule_index'):
				continue
			elif isinstance(self,Atom) and a=='bonded':
				if len(selfDict[a]) == len(otherDict[a]):
					continue
				else:
					return False

			#Checking _Physical attributes.
			elif isinstance(selfDict[a],_Physical):
				if selfDict[a].equals(otherDict[a]):
					continue
				else:
					return False

				#If an attribute is a list or tuple, go through it element-by-element,
				#assuming the order is the same between other and self, and if
				#any of the lists's elements are _Physical instances, compare them
				#with .equals()

				#Exception: If self is an Atom, do not compare the bonded lists.
				#This would result in an infinite loop.
			elif isinstance(selfDict[a],list) or isinstance(selfDict[a],
					tuple):
				if len(selfDict[a])!= len(otherDict[a]):
					return False
				#print "Self "+ `type(selfDict[a])`+ " is "+`selfDict[a]`
				#print "Other " + `type(otherDict[a])`+ " is "+`otherDict[a]`
				for b in range(len(selfDict[a])):
					if selfDict[a][b] == otherDict[a][b]:
						continue
					elif isinstance(selfDict[a][b],_Physical):
						if selfDict[a][b].equals(otherDict[a][b]):
							continue
						else:
							#print "False 1"
							return False
					else:
						#print "False 2"
						return False
			else:
				#print "False 3"
				return False
		#print "Passed All Check Blocks. Returning True."
		return True

	def __repr__(self):
		text = self.__dict__
		if "bonded" in text:
			del text["bonded"]
		return object.__repr__(self) +" with attributes:\n"+str(text)

class Atom(_Physical):
	def __init__(self, element, x, y, z, vx=0, vy=0, vz=0, fx=0, fy=0, fz=0, index=None, type=None, molecule_index=1, bonded=[], type_index=None):
		self.element = element
		self.x = x
		self.y = y
		self.z = z
		self.vx = vx
		self.vy = vy
		self.vz = vz
		self.fx = fx
		self.fy = fy
		self.fz = fz
		self.index = index
		self.molecule_index = molecule_index
		self.bonded=bonded
		self.type_index=type_index
		# When parameterized by OPLS, 'type' dict contains: {'bond_count': 3, 'index': 588, 'notes': '1,10-Phenanthroline C2', 'vdw_r': 3.55, 'element': 6, 'vdw_e': 0.07, 'charge': 0.392, 'mass': 12.011, 'index2': 48, 'element_name': 'CA'}
		# Object May also contain lammps_type, mass, and charge
		self.type = type

	def translate(self, v):
		self.x+=v[0]; self.y+=v[1]; self.z+=v[2]

	def to_string(self, verbose=False):
		if self.index == None:
			raise ValueError("Atom cannot have index 'None'.")
		text = '%s, (%3.3f, %3.3f, %3.3f), index: %d' % (self.element, self.x, self.y, self.z, self.index)

		if self.type_index:
			text += ', type_index: %d' % (self.type_index)

		text += ', molecule_index: %d' % (self.molecule_index)

		if self.type and verbose:
			text += str(self.type)

		text += '\n'

		return text

	def flatten(self):
		return [self.x,self.y,self.z]

	def set_position(self,pos):
		self.x = pos[0]
		self.y = pos[1]
		self.z = pos[2]

	def __str__(self):
		return self.to_string()

class Bond(_Physical):
	def __init__(self, a, b, type=None, r=None):
		self.atoms = (a,b)
		self.type = type
		self.r = r

class Angle(_Physical):
	def __init__(self, a, b, c, type=None, theta=None):
		self.atoms = (a,b,c)
		self.type = type
		self.theta = theta

	def __repr__(self):
		text = self.__dict__
		return object.__repr__(self) +" with attributes:\n"+str(text)

class Dihedral(_Physical):
	def __init__(self, a, b, c, d, type=None, theta=None):
		self.atoms = (a,b,c,d)
		self.type = type
		self.theta = theta

class Job(object): #a job on the queue
	def __init__(self, name):
		self.name = name
	def wait(self):
		while True:
			jlist = subprocess.Popen('jlist', shell=True, stdout=subprocess.PIPE).communicate()[0]
			if (' '+self.name+' ') in jlist:
				time.sleep(60)
			else: break

# A generic class to hold dft data
class DFT_out(object):
	def __init__(self, name, dft='g09'):
		self.name = name
		self.dft = dft.lower()

		# Initialize everything as None
		self.route = None
		self.frames = None
		self.atoms = None
		self.gradients=None
		self.energies = None
		self.energy = None
		self.charges_MULLIKEN = None
		self.charges_LOEWDIN = None
		self.charges_CHELPG = None
		self.charges = None
		self.convergence = None
		self.converged = None
		self.time = None
		self.bandgaps = None
		self.bandgap = None
		self.finished = None
		self.warnings = None

# A generic class to hold simulation data, aprticularly lammps trajectory files
class sim_out(object):
	def __init__(self, name, program='lammps'):
		self.name = name
		self.program = program.lower()

		# Initialize everything as None
		self.frames = None
		self.atoms = None
		self.timesteps = None
		self.final_timestep = None
		self.atom_counts = None
		self.atom_count = None
		self.box_bounds_list = None
		self.box_bounds = None

def get_bonds(atoms):
	bonds = []
	for a in atoms:
		a.bonded = []
	for i,a in enumerate(atoms):
		for b in atoms[i+1:]:
			dd = dist_squared(a,b)
			if (a.element not in [1,'H'] and b.element not in [1,'H'] and dd<2**2) or (dd < 1.2**2 and (a.element in [1,'H'])!=(b.element in [1,'H']) ) or (dd < 3.1**2 and (a.element in ['Pb',82] or b.element in ['Pb',82]) ):
				bonds.append( Bond(a,b,r=dd**0.5) )
				if a not in b.bonded: b.bonded.append(a)
				if b not in a.bonded: a.bonded.append(b)
	return bonds

def get_angles_and_dihedrals(atoms):
	angles = []
	for center in atoms:
		if len(center.bonded)<2: continue
		for i,a in enumerate(center.bonded):
			for b in center.bonded[i+1:]:
				A = math.sqrt((center.z-b.z)**2+(center.x-b.x)**2+(center.y-b.y)**2)
				N = math.sqrt((a.z-b.z)**2+(a.x-b.x)**2+(a.y-b.y)**2)
				B = math.sqrt((center.z-a.z)**2+(center.x-a.x)**2+(center.y-a.y)**2)
				try:
					theta = 180/math.pi*math.acos((A**2+B**2-N**2)/(2*A*B))
				except: theta = 0.0
				angles.append( Angle(a,center,b,theta=theta) )
	#Updated to provide deterministic dihedral order with the same time complexity
	dihedral_list = []
	dihedral_set = {}
	for angle in angles:
		for a in angle.atoms[0].bonded:
			if a is angle.atoms[1]: continue
			dihedral = (a,) + angle.atoms
			if tuple(reversed(dihedral)) not in dihedral_set:
				dihedral_set[dihedral] = True
				dihedral_list.append(dihedral)

		for b in angle.atoms[2].bonded:
			if b is angle.atoms[1]: continue
			dihedral = angle.atoms + (b,)
			if tuple(reversed(dihedral)) not in dihedral_set:
				dihedral_set[dihedral] = True
				dihedral_list.append(dihedral)
	dihedral_list = _Uniquify(dihedral_list)
	dihedrals = [Dihedral(*d) for d in dihedral_list]

	return angles, dihedrals

def frange(low, high, step):
	warn(
        "similar functionality already exists in numpy. Please switch to np.arange.",
        DeprecationWarning
    )
	while low < high:
		yield low
		low += step

def quat_to_mat(q): #quat = [w i j k]
	warn(
        "this functionality may be possible in numpy.",
        DeprecationWarning
    )
	d,b,c,a = q
	return [ [a**2+b**2-c**2-d**2, 2*b*c-2*a*d, 2*b*d+2*a*c],
	[2*b*c+2*a*d, a**2-b**2+c**2-d**2, 2*c*d-2*a*b],
	[2*b*d-2*a*c, 2*c*d-2*a*b, a**2-b**2-c**2+d**2]
	]

def matvec(m,v):
	warn(
        "this functionality may be possible in numpy.",
        DeprecationWarning
    )
	return (m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2], m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2], m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2])

def matmat(a,b):
	warn(
        "this functionality may be possible in numpy.",
        DeprecationWarning
    )
	product = [[0.]*3, [0.]*3, [0.]*3]
	for y in range(3):
		for x in range(3):
			for k in range(3):
				product[y][x] += a[y][k]*b[k][x]
	return product


def rand_rotation(): #http://tog.acm.org/resources/GraphicsGems/, Ed III
	x = (random.random(), random.random(), random.random())
	theta = x[0] * 2*math.pi
	phi   = x[1] * 2*math.pi
	z = x[2] * 2
	#Compute a vector V used for distributing points over the sphere via the reflection I - V Transpose(V).  This formulation of V will guarantee that if x[1] and x[2] are uniformly distributed, the reflected points will be uniform on the sphere.  Note that V has length sqrt(2) to eliminate the 2 in the Householder matrix.
	r = math.sqrt(z)
	Vx = math.sin( phi ) * r
	Vy = math.cos( phi ) * r
	Vz = math.sqrt( 2.0 - z )
	#Compute the row vector S = Transpose(V) * R, where R is a simple rotation by theta about the z-axis.  No need to compute Sz since it's just Vz.
	st = math.sin( theta )
	ct = math.cos( theta )
	Sx = Vx * ct - Vy * st
	Sy = Vx * st + Vy * ct

	#Construct the rotation matrix  ( V Transpose(V) - I ) R, which is equivalent to V S - R.
	M = [ [0.,0.,0.], [0.,0.,0.], [0.,0.,0.] ]

	M[0][0] = Vx * Sx - ct
	M[0][1] = Vx * Sy - st
	M[0][2] = Vx * Vz

	M[1][0] = Vy * Sx + st
	M[1][1] = Vy * Sy - ct
	M[1][2] = Vy * Vz

	M[2][0] = Vz * Sx
	M[2][1] = Vz * Sy
	M[2][2] = 1.0 - z	# This equals Vz * Vz - 1.0

	return M

# Construct rotation matrix around x, using t
def rot_x(t):
	M = [ [1., 0., 0.], [0., math.cos(t), -math.sin(t)], [0., math.sin(t), math.cos(t)] ]
	return M

# Construct rotation matrix around y, using t
def rot_y(t):
	M = [ [math.cos(t), 0., math.sin(t)], [0., 1., 0.], [-math.sin(t), 0., math.cos(t)] ]
	return M

# Construct rotation matrix around z, using t
def rot_z(t):
	M = [ [math.cos(t), -math.sin(t), 0.], [math.sin(t), math.cos(t), 0.], [0., 0., 1.] ]
	return M

# Construct general rotation matrix using yaw, pitch, and roll (alpha, beta, gamma)
# Performs extrinsic rotation whose Euler angles are alpha, beta, gamma about axes z, y, x
def rot_xyz(alpha, beta, gamma):
	# Extrinsic definition
	M = matmat(matmat(rot_z(gamma), rot_y(beta)), rot_x(alpha))
	# Intrinsic definition
	#M = matmat(matmat(rot_z(alpha), rot_y(beta)), rot_x(gamma))
	return M

class Molecule(_Physical):
	def __init__(self, atoms_or_filename_or_all, bonds=None, angles=None, dihedrals=None, parameter_file=opls_path, extra_parameters={}, test_charges=True, allow_errors=False): #set atoms, bonds, etc, or assume 'atoms' contains all those things if only one parameter is passed in
		if type(atoms_or_filename_or_all)==type('string'):
			self.filename = atoms_or_filename_or_all
			atoms, bonds, angles, dihedrals = files.read_cml(self.filename, parameter_file=parameter_file, extra_parameters=extra_parameters, test_charges=test_charges, allow_errors=allow_errors)
		elif bonds:
			atoms, bonds, angles, dihedrals = atoms_or_filename_or_all
		else:
			atoms, bonds, angles, dihedrals = atoms_or_filename_or_all, [], [], []
		self.atoms = atoms
		self.bonds = bonds
		self.angles = angles
		self.dihedrals = dihedrals
	def rotate(self, m):
		for a in self.atoms:
			a.x, a.y, a.z = np.dot(np.asarray(m),np.array([a.x,a.y,a.z]))

	def randRotateInPlace(self):
		"""Randomly rotates the molecule around its center of mass"""
		center = self.CalcCenter()
		self.SetCenter([0.0,0.0,0.0])
		self.rand_rotate()
		self.translate(center)

	def translate(self, v):
		for a in self.atoms:
			a.x+=v[0]; a.y+=v[1]; a.z += v[2]
	def rand_rotate(self):
		rand_m = rand_rotation()
		self.rotate(rand_m)

	def get_com(self, skip_H=True):
		if skip_H: n = float(len([a for a in self.atoms if a.element != "H"]))
		else: n = float(len(self.atoms))
		if skip_H:
			x = sum([a.x for a in self.atoms if a.element != "H"])/n
			y = sum([a.y for a in self.atoms if a.element != "H"])/n
			z = sum([a.z for a in self.atoms if a.element != "H"])/n
		else:
			x = sum([a.x for a in self.atoms])/n
			y = sum([a.y for a in self.atoms])/n
			z = sum([a.z for a in self.atoms])/n
		return (x,y,z)

	# Print all atoms
	def to_string(self):
		text = ''
		for atom in self.atoms:
			text += atom.to_string()
		return text

	def flatten(self):
		return np.array([[a.x, a.y, a.z] for a in self.atoms]).flatten()

	def set_positions(self, positions, new_atom_list=False):
		positions = positions.flatten().reshape((-1,3))
		if len(positions) != len(self.atoms) and not new_atom_list:
			raise Exception("position list does not hold the same number of atoms as does this molecule. Consider raising new_atom_list flag in set_positions.")
		if new_atom_list: self.atoms = [Atom("",p[0],p[1],p[2]) for p in positions]
		else:
			for a,b in zip(self.atoms, positions):
				a.x, a.y, a.z = b[0], b[1], b[2]

	def CalcCenter(self):
		"""Returns the center of mass of the molecule in [x,y,z] format."""
		xList = []
		yList = []
		zList = []
		totalMass = 0.0
		for a in self.atoms:
			xList.append(a.x*a.type.mass)
			yList.append(a.y*a.type.mass)
			zList.append(a.z*a.type.mass)
			totalMass += a.type.mass

		return [sum(xList)/totalMass,sum(yList)/totalMass,sum(zList)/totalMass]

	def SetCenter(self,xyz=[0.0,0.0,0.0]):
		"""
		Sets the center of mass of the molecule to the given xyz position,
		in [x,y,z] format. If no xyz is passed, resets the center of mass to the
		origin.
		"""
		center = self.CalcCenter()
		self.translate([xyz[0]-center[0],xyz[1]-center[1],xyz[2]-center[2]])

	# When printing molecule, print all atoms
	def __str__(self):
		return self.to_string()



class System(_Physical):
	# Initialize all system variables. Convert to float as necessary to prepare for any needed math
	def __init__(self, name=None, box_size=(10.0,10.0,10.0), box_angles=(90.0,90.0,90.0), periodic=False):
		self.atoms, self.bonds, self.angles, self.dihedrals = [], [], [], []
		self.box_size = [float(param) for param in box_size] # (a, b, c)
		self.box_angles = [float(param) for param in box_angles] # (alpha, beta, gamma)
		self.name = name
		self.molecules = []
		self.periodic = periodic

		# If the system is not a monoclinic box, set lammps triclinic parameters
		if abs(box_angles[0]-90) > 0.001 or abs(box_angles[1]-90) > 0.001 or abs(box_angles[2]-90) > 0.001:
			self.setTriclinicBox(self.periodic,self.box_size,self.box_angles)

		# If system is a monoclinic box and set default lammps triclinic parameters
		# Assumes center of box is the origin
		else:
			self.xlo = -self.box_size[0]/2.0
			self.ylo = -self.box_size[1]/2.0
			self.zlo = -self.box_size[2]/2.0
			self.xhi = self.box_size[0]/2.0
			self.yhi = self.box_size[1]/2.0
			self.zhi = self.box_size[2]/2.0
			self.xy = 0.0
			self.xz = 0.0
			self.yz = 0.0

	# Establish lammps triclinic box boundary conditions for this system
	def setTriclinicBox(self, periodic, box_size, box_angles):
		self.box_size[0] = box_size[0]
		self.box_size[1] = box_size[1]
		self.box_size[2] = box_size[2]
		self.box_angles[0] = box_angles[0]
		self.box_angles[1] = box_angles[1]
		self.box_angles[2] = box_angles[2]

		a = box_size[0]
		b = box_size[1]
		c = box_size[2]
		alpha = box_angles[0]
		beta = box_angles[1]
		gamma = box_angles[2]

		# For lammmps, trigonal vectors established by using xy xz yz
		# A = (xhi-xlo,0,0);
		# B = (xy,yhi-ylo,0);
		# C = (xz,yz,zhi-zlo)
		self.xlo = 0.0
		self.ylo = 0.0
		self.zlo = 0.0

		# Formula for converting (a,b,c,alpha,beta,gamma) to (lx,ly,lz,xy,xz,yz)
		# taken from online lammps help
		self.xhi = a
		self.xy = b*math.cos(math.radians(gamma))
		self.xz = c*math.cos(math.radians(beta))
		self.yhi = math.sqrt(b**2 - self.xy**2)
		self.yz = (b*c*math.cos(math.radians(alpha)) - self.xy * self.xz)/ self.yhi
		self.zhi = math.sqrt(c**2 - self.xz**2 - self.yz**2)

	def add(self, molecule, x=0.0, y=0.0, z=0.0, scale_x=1, scale_y=1, scale_z=1):
		atom_offset = len(self.atoms)
		for a in molecule.atoms:
			new_atom = copy.copy(a)
			new_atom.index=a.index+atom_offset
			new_atom.x*=scale_x; new_atom.y*=scale_y; new_atom.z*=scale_z
			new_atom.x+=x; new_atom.y+=y; new_atom.z+=z
			new_atom.molecule_index=len(self.molecules)+1
			self.atoms.append( new_atom )
		for t in molecule.bonds:
			self.bonds.append( Bond(*[self.atoms[a.index+atom_offset-1] for a in t.atoms], type=t.type) )
		for t in molecule.angles:
			self.angles.append( Angle(*[self.atoms[a.index+atom_offset-1] for a in t.atoms], type=t.type) )
		for t in molecule.dihedrals:
			self.dihedrals.append( Dihedral(*[self.atoms[a.index+atom_offset-1] for a in t.atoms], type=t.type) )
		new_molecule = copy.copy(molecule)
		new_molecule.atoms = self.atoms[-len(molecule.atoms):]
		new_molecule.bonds = self.bonds[-len(molecule.bonds):]
		new_molecule.angles = self.angles[-len(molecule.angles):]
		new_molecule.dihedrals = self.dihedrals[-len(molecule.dihedrals):]
		self.molecules.append( new_molecule )

	def Remove(self,target):
		"""
		If target is a molecule, removes all atoms, bonds angles and dihedrals of
		the passed molecule from the system. Raises a ValueError if not all aspects
		of molecule are found in the system.

		If target is an Atom, the atom is removed from the system, and any bonds,
		angles, and dihedrals which contain the atom are also removed from the system.
		Raises a ValueError if the Atom is not found in the system.

		Precondition: molecule is a valid Utils.Molecule instance or a valid
		Utils.Atom instance.
		"""
		#Make sure all atoms in molecule are in system.
		if isinstance(target,Molecule):
			if len(target.atoms)>0:
				for a in target.atoms:
					atomCheck = False
					for b in range(len(self.atoms)):
						#Check from the back of the atomlist, because the atom
						#to be removed was also likely the last one added.
						if self.atoms[len(self.atoms)-b-1].equals(a):
							del self.atoms[len(self.atoms)-b-1]
							atomCheck = True
							break
					if not atomCheck:
						raise ValueError("_Physical instance "+`a`+" wasn't found"+
										 "in the given system.")

			#Repeat above for bonds, angles, dihedrals
			if len(target.bonds)>0:
				for a in target.bonds:
					bondCheck = False
					for b in range(len(self.bonds)):
						if self.bonds[len(self.bonds)-b-1].equals(a):
							bondCheck = True
							del self.bonds[len(self.bonds)-b-1]
							break
					if not bondCheck:
						raise ValueError("_Physical instance "+`a`+" wasn't found"+
										 "in the given system.")

			if len(target.angles)>0:
				for a in target.angles:
					angleCheck = False
					for b in range(len(self.angles)):
						if self.angles[len(self.angles)-b-1].equals(a):
							angleCheck = True
							del self.angles[len(self.angles)-b-1]
							break
					if not angleCheck:
						raise ValueError("_Physical instance "+`a`+" wasn't found"+
										 "in the given system.")

			if len(target.dihedrals)>0:
				for a in target.dihedrals:
					dihedralCheck = False
					for b in range(len(self.dihedrals)):
						if self.dihedrals[len(self.dihedrals)-b-1].equals(a):
							dihedralCheck = True
							del self.dihedrals[len(self.dihedrals)-b-1]
							break
					if not dihedralCheck:
						raise ValueError("_Physical instance "+`a`+" wasn't found"+
										 "in the given system.")

			for a in range(len(self.molecules)):
				if self.molecules[a].equals(target):
					del self.molecules[a]
					break

		elif isinstance(target,Atom):
			for a in range(len(self.atoms)):
				#Check from the back of the atomlist, because the atom
				#to be removed was also likely the last one added.
				atomCheck = False
				if target.equals(self.atoms[-a-1]):
					del self.atoms[-a-1]
					atomCheck = True
					break
			if not atomCheck:
				raise ValueError("_Physical instance "+`target`+" wasn't found"+
									 "in the given system.")

			delList= []
			newList= []
			for a in range(len(self.bonds)):
				for b in self.bonds[a].atoms:

					#If a bond has the target atom, mark it
					if b.equals(target):
						delList.append(a)

			#Delete all of the bonds that were marked, from the end of the list
			#to the front so as to preserve order
			for a in range(len(self.bonds)):
				if not a in delList:
					newList.append(self.bonds[a])
			self.bonds=newList

			delList= []
			newList= []
			for a in range(len(self.angles)):
				for b in self.angles[a].atoms:
					if b.equals(target):
						delList.append(a)
			for a in range(len(self.angles)):
				if not a in delList:
					newList.append(self.angles[a])
			self.angles=newList

			delList= []
			newList= []
			for a in range(len(self.dihedrals)):
				for b in self.dihedrals[a].atoms:
					if b.equals(target):
						delList.append(a)
			for a in range(len(self.dihedrals)):
				if not a in delList:
					newList.append(self.dihedrals[a])
			self.dihedrals=newList

	def Contains(self,molecule):
		"""
		Returns a boolean, True if the molecule passed as an argument is contained
		exactly within this system. This implies that all atoms,bonds,angles and
		dihedrals in molecule are in the system. Matching is checked using the
		_Physical.equals() method. Returns false if it is not the case that
		all properties of the molecule are found exactly as given in the current
		system.
		"""
		#Make sure all atoms in molecule are in system.
		if len(molecule.atoms)>0:
			for a in molecule.atoms:
				atomCheck = False
				for b in range(len(self.atoms)):
					if self.atoms[len(self.atoms)-b-1].equals(a):
						atomCheck = True
						break
				if not atomCheck:
					return False

		#Repeat above for bonds, angles, dihedrals
		if len(molecule.bonds)>0:
			for a in molecule.bonds:
				bondCheck = False
				for b in range(len(self.bonds)):
					if self.bonds[len(self.bonds)-b-1].equals(a):
						bondCheck = True
						break
				if not bondCheck:
					return False

		if len(molecule.angles)>0:
			for a in molecule.angles:
				angleCheck = False
				for b in range(len(self.angles)):
					if self.angles[len(self.angles)-b-1].equals(a):
						angleCheck = True
						break
				if not angleCheck:
					return False

		if len(molecule.dihedrals)>0:
			for a in molecule.dihedrals:
				dihedralCheck = False
				for b in range(len(self.dihedrals)):
					if self.dihedrals[len(self.dihedrals)-b-1].equals(a):
						dihedralCheck = True
						break
				if not dihedralCheck:
					return False
		#If python gets here, all aspects of molecule are in system.
		return True


	# Print all atoms
	def to_string(self):
		text = ''
		for atom in self.atoms:
			text += atom.to_string()
		return text


	def __str__(self):
		return self.to_string()


	def __repr__(self):
		text = ''
		for atom in self.atoms:
			text += atom.to_string()
		return text


def _Uniquify(givenList,idfun=None):
	"""
	Returns the given list with order preserved, with duplicates past the first
	entry removed.
	Credits to: https://www.peterbe.com/plog/uniqifiers-benchmark
	"""
	# Order preserving
	return list(_f10(givenList, idfun))

def _f10(seq, idfun=None):
	"""Helper function to _Uniquify()"""
	seen = set()
	if idfun is None:
		for x in seq:
			if x in seen:
				continue
			seen.add(x)
			yield x
	else:
		for x in seq:
			x = idfun(x)
			if x in seen:
				continue
			seen.add(x)
			yield x

def dist_squared(a,b):
	return (a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2

def dist(a,b):
	return dist_squared(a,b)**0.5

def angle_size(a,center,b):
	A = math.sqrt((center.z-b.z)**2+(center.x-b.x)**2+(center.y-b.y)**2)
	N = math.sqrt((a.z-b.z)**2+(a.x-b.x)**2+(a.y-b.y)**2)
	B = math.sqrt((center.z-a.z)**2+(center.x-a.x)**2+(center.y-a.y)**2)
	return 180/math.pi*math.acos((A**2+B**2-N**2)/(2*A*B))

def dihedral_angle(a,b,c,d):
    """Praxeolitic formula
    1 sqrt, 1 cross product
    http://stackoverflow.com/a/34245697"""
    from numpy import array, dot, degrees, arctan2, cross, cos
    from numpy.linalg import norm
    p0 = array([a.x,a.y,a.z])
    p1 = array([b.x,b.y,b.z])
    p2 = array([c.x,c.y,c.z])
    p3 = array([d.x,d.y,d.z])

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - dot(b0, b1)*b1
    w = b2 - dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = dot(v, w)
    y = dot(cross(b1, v), w)

    phi = arctan2(y, x)
    return phi, cos(phi), cos(2*phi), cos(3*phi), cos(4*phi)

# A test procrustes code to remove possibility of reflection
def orthogonal_procrustes(A, ref_matrix, reflection=False):
	# Adaptation of scipy.linalg.orthogonal_procrustes -> https://github.com/scipy/scipy/blob/v0.16.0/scipy/linalg/_procrustes.py#L14
	# Info here: http://compgroups.net/comp.soft-sys.matlab/procrustes-analysis-without-reflection/896635
	# goal is to find unitary matrix R with det(R) > 0 such that ||A*R - ref_matrix||^2 is minimized
	from scipy.linalg.decomp_svd import svd # Singular Value Decomposition, factors matrices
	from scipy.linalg import det
	import numpy as np

	A = np.asarray_chkfinite(A)
	ref_matrix = np.asarray_chkfinite(ref_matrix)

	if A.ndim != 2:
		raise ValueError('expected ndim to be 2, but observed %s' % A.ndim)
	if A.shape != ref_matrix.shape:
		raise ValueError('the shapes of A and ref_matrix differ (%s vs %s)' % (A.shape, ref_matrix.shape))


	u, w, vt = svd(ref_matrix.T.dot(A).T)

	# Goal: minimize ||A*R - ref||^2, switch to trace
	# trace((A*R-ref).T*(A*R-ref)), now we distribute
	# trace(R'*A'*A*R) + trace(ref.T*ref) - trace((A*R).T*ref) - trace(ref.T*(A*R)), trace doesn't care about order, so re-order
	# trace(R*R.T*A.T*A) + trace(ref.T*ref) - trace(R.T*A.T*ref) - trace(ref.T*A*R), simplify
	# trace(A.T*A) + trace(ref.T*ref) - 2*trace(ref.T*A*R)
	# Thus, to minimize we want to maximize trace(ref.T * A * R)

	# u*w*v.T = (ref.T*A).T
	# ref.T * A = w * u.T * v
	# trace(ref.T * A * R) = trace (w * u.T * v * R)
	# differences minimized when trace(ref.T * A * R) is maximized, thus when trace(u.T * v * R) is maximized
	# This occurs when u.T * v * R = I (as u, v and R are all unitary matrices so max is 1)
	# R is a rotation matrix so R.T = R^-1
	# u.T * v * I = R^-1 = R.T
	# R = u * v.T
	# Thus, R = u.dot(vt)

	R = u.dot(vt) # Get the rotation matrix, including reflections
	if not reflection and det(R) < 0: # If we don't want reflection
		# To remove reflection, we change the sign of the rightmost column of u (or v) and the scalar associated
		# with that column
		u[:,-1] *= -1
		w[-1] *= -1
		R = u.dot(vt)

	scale = w.sum() # Get the scaled difference

	return R,scale

# Procrustes works by geting an orthogonal frame to map frames[1:] to be as similar to frames[0] as possible
# This implements the orthagonal procrustes with translation and no reflection (Partial Procrustes)
def procrustes(frames, count_atoms=None, append_in_loop=True):
	if not count_atoms: count_atoms = range(len(frames[0]))
	for s in frames:
		center_x = sum([a.x for i,a in enumerate(s) if i in count_atoms])/len(count_atoms)
		center_y = sum([a.y for i,a in enumerate(s) if i in count_atoms])/len(count_atoms)
		center_z = sum([a.z for i,a in enumerate(s) if i in count_atoms])/len(count_atoms)
		for a in s:
			a.x -= center_x
			a.y -= center_y
			a.z -= center_z
	# rotate all frames to be as similar to their neighbors as possible
	from scipy.linalg import det
	from numpy import dot

	full_rotation = []

	for i in range(1,len(frames)): #rotate all frames to optimal alignment
		# only count spring-held atoms for finding alignment
		# orthogonal_procrustes maps count_atoms_1 onto count_atoms_2
		count_atoms_1 = [(a.x,a.y,a.z) for j,a in enumerate(frames[i]) if j in count_atoms]
		count_atoms_2 = [(a.x,a.y,a.z) for j,a in enumerate(frames[i-1]) if j in count_atoms]
		rotation = orthogonal_procrustes(count_atoms_1,count_atoms_2)[0]


		if det(rotation) < 0:
			raise Exception('Procrustes returned reflection matrix')
		# rotate all atoms into alignment
		for a in frames[i]:
			a.x,a.y,a.z = dot((a.x,a.y,a.z), rotation)
			if hasattr(a,'fx'): a.fx,a.fy,a.fz = dot((a.fx,a.fy,a.fz), rotation)
			if append_in_loop:
				full_rotation.append(rotation)
		if not append_in_loop:
			full_rotation.append(rotation)

	return full_rotation

def interpolate(atoms1, atoms2, N): #interpolate N steps between two sets of coordinates
	frames = [[] for i in range(N)]
	for a,b in zip(atoms1,atoms2):
		dx,dy,dz = b.x-a.x, b.y-a.y, b.z-a.z
		for i in range(N):
			frac = 1.0*(i+1)/(N+1)
			frames[i].append( Atom(a.element, a.x+dx*frac, a.y+dy*frac, a.z+dz*frac) )
	return frames

def motion_per_frame(frames):
	per_state_avg = [0.0 for s in frames]
	for atom_list in zip(*frames):
		for i in range(1,len(atom_list)):
			a = atom_list[i-1]
			b = atom_list[i]
			per_state_avg[i] += dist(a,b)
	motion = []
	for x in per_state_avg: motion.append(x/len(frames[0]))
	return motion

# Alternative to Procrustes.  This code takes 3 values in a list as ids to recenter frames:
# ids[0] - This is an atom that will be positioned at the origin after translating the frame
# ids[1] - This is an atom that will lie on the positive x-axis after two rotations of the frame
# ids[2] - This is an atom that will lie on the xy plane in the positive y direction after rotation of the frame
def center_frames(frames,ids,X_TOL=0.1,XY_TOL=0.1,Z_TOL=0.1,THETA_STEP=0.005,TRANSLATE=[0,0,0]):
	from math import sin, cos, pi

	def get_pnt(a): return [a.x,a.y,a.z]
	def trans(a,t):
		try:
			a.x += t.x
			a.y += t.y
			a.z += t.z
		except:
			a.x += t[0]
			a.y += t[1]
			a.z += t[2]
	def rot_yz(a,t):
		x = a.x
		y = a.y*cos(t)-a.z*sin(t)
		z = a.y*sin(t)+a.z*cos(t)
		a.x = x
		a.y = y
		a.z = z

		try:
			fx = a.fx
			fy = a.fy*cos(t)-a.fz*sin(t)
			fz = a.fy*sin(t)+a.fz*cos(t)
			a.fx = fx
			a.fy = fy
			a.fz = fz
		except: pass

	def rot_xy(a,t):
		x = a.x*cos(t)-a.y*sin(t)
		y = a.x*sin(t)+a.y*cos(t)
		z = a.z
		a.x = x
		a.y = y
		a.z = z

		try:
			fx = a.fx*cos(t)-a.fy*sin(t)
			fy = a.fx*sin(t)+a.fy*cos(t)
			fz = a.fz
			a.fx = fx
			a.fy = fy
			a.fz = fz
		except: pass

	origin = ids[0]
	xaxis = ids[1]
	sqr = ids[2]

	# If we only have one frame, put it in a list
	chk = False
	if type(frames[0]) != list:
		frames = [frames]
		chk = True

	# Loop through frames
	for f in frames:
		# Find the first translation to make the desired point the origin
		trans_1 = get_pnt(f[origin])
		for i in range(len(trans_1)): trans_1[i] *= -1

		# Translate everything
		for a in f: trans(a,trans_1)

		# Find the desired x-axis' rotation to place it on the xy plane
		theta = 0
		pnt = f[xaxis]
		while 1:
			chk = pnt.y*sin(theta) + pnt.z*cos(theta)
			if abs(chk) < Z_TOL: break
			theta += THETA_STEP
			if theta > 2*pi:
				print("Cannot place atom of index %d on the xy plane" % xaxis)
				print("Try decreasing THETA_STEP below %lg..." % THETA_STEP)
				sys.exit()

		# Rotate everything
		for a in f: rot_yz(a,theta)

		# Now find the angle that we rotate around the z axis to get the +x-axis aligned
		theta = 0
		pnt = f[xaxis]
		while 1:
			chk_x = pnt.x*cos(theta) - pnt.y*sin(theta)
			chk = pnt.x*sin(theta) + pnt.y*cos(theta)
			if abs(chk) < X_TOL and chk_x > 0: break
			theta += THETA_STEP
			if theta > 2*pi:
				print("Cannot place atom of index %d on the x axis" % xaxis)
				print("Try decreasing THETA_STEP below %lg..." % THETA_STEP)
				sys.exit()

		# Rotate everything
		for a in f: rot_xy(a,theta)

		# Now find the angle that we rotate around the x axis such that our last vector lies on the x(+y) plane
		theta = 0
		pnt = f[sqr]
		while 1:
			chk_y = pnt.y*cos(theta)-pnt.z*sin(theta)
			chk = pnt.y*sin(theta) + pnt.z*cos(theta)
			if abs(chk) < XY_TOL and chk_y > 0: break
			theta += THETA_STEP
			if theta > 2*pi:
				print("Cannot place atom of index %d on the x(+y) plane" % sqr)
				print("Try decreasing THETA_STEP below %lg..." % THETA_STEP)
				sys.exit()

		# Rotate everything
		for a in f: rot_yz(a,theta)

		# Re-translate the whole system
		for a in f: trans(a,TRANSLATE)

	from math import isnan
	for f in frames:
		for a in f:
			if isnan(a.x) or isnan(a.y) or isnan(a.z):
				print("Center frames has led to NaN...")
				sys.exit()

	if chk: frames = frames[0]

def pretty_xyz(name,R_MAX=1,F_MAX=50,PROCRUSTES=False,outName=None,write_xyz=False,verbose=False):
	#----------
	# name = Name of xyz file to read in
	# R_MAX = maximum motion per frame
	# F_MAX = maximum number of frames allowed
	# PROCRUSTES = Center frames or not
	# outName = If you wish to output to file, give a name. Defalt name is 'pretty_xyz'
	# write_xyz = Write to file. Default False
	# Verbose = Outputing what pretty_xyz is doing as it goes
	#----------

	from copy import deepcopy

	# Get data as either frames or a file
	if type(name)==type(''): frames = files.read_xyz(name)
	elif type(name)==type([]): frames = name
	else:
		print "Error - Invalid name input.  Should be either the name of an xyz file or a list.", sys.exc_info()[0]
		exit()

	# Loop till we're below R_MAX
	while 1:
		# Find largest motion_per_frame
		if PROCRUSTES: procrustes(frames)
		tmp = motion_per_frame(frames)
		i = tmp.index(max(tmp))

		# Check if we're done
		r2 = max(tmp)
		if r2 < R_MAX: break

		if len(frames) > F_MAX:
			print "-------------------------------------------------------"
			print tmp
			print "-------------------------------------------------------"
			print "\n\nError - Could not lower motion below %lg in %d frames." % (R_MAX,F_MAX), sys.exc_info()[0]
			exit()
		else:
			if verbose: print "Currently Frames = %d\tr2 = %lg" % (len(frames),r2)

		# Now, split the list, interpolate, and regenerate
		if i > 0 and i < len(frames) - 1:
			f_low = deepcopy(frames[:i])
			f_high = deepcopy(frames[i+1:])
			f_mid = interpolate(frames[i-1],frames[i+1],3)
			frames = f_low + f_mid + f_high
		elif i == 0:
			f_low = deepcopy(frames[i])
			f_mid = interpolate(frames[i],frames[i+1],3)
			f_high = deepcopy(frames[i+1:])
			frames = [f_low] + f_mid + f_high
		else:
			f_low = deepcopy(frames[:i])
			f_mid = interpolate(frames[i-1],frames[i],3)
			f_high = deepcopy(frames[i])
			frames = f_low + f_mid + [f_high]

		if verbose: print "\tInterpolated %d,%d ... %lg" % (i-1,i+1,max(motion_per_frame(frames)))

	if PROCRUSTES: procrustes(frames)
	if verbose: print("\tThere are now a total of %d frames" % len(frames))

	if write_xyz: files.write_xyz(frames,'pretty_xyz' if outName==None else outName)
	else: return frames

# A function to format a string's colour
def color_set(s,c): return constants.COLOR[c] + str(s) + constants.COLOR['ENDC']
colour_set = color_set

def strip_color(s):
	for c in constants.COLOUR:
		col = constants.COLOUR[c]
		while col in s:
			s = s.replace(col,'')
	return s
strip_colour = strip_color

def opls_options(molecule, parameter_file=opls_path):
	elements, atom_types, bond_types, angle_types, dihedral_types = files.read_opls_parameters(parameter_file)

	elements_by_structure_indices = dict( [ (t.index2, elem_i2s(t.element) ) for t in atom_types ] )
	elements_by_structure_indices[0] = 'X'

	def add_to_list(dic,key,value):
		if key in dic:
			dic[key].append(value)
		else:
			dic[key] = [value]

	dihedral_types_by_element={}
	for d in dihedral_types:
		structure_indices = d.index2s
		elements = [ elements_by_structure_indices[i] for i in structure_indices]
		add_to_list(dihedral_types_by_element, tuple(elements), d)
		add_to_list(dihedral_types_by_element, tuple(reversed(elements)), d)

	atoms, bonds, angles, dihedrals = files.read_cml(molecule, parameter_file=None)

	for a in atoms:
		a.index2_options = []

	for d in dihedrals:
		elements = tuple([ a.element for a in d.atoms ])
		options = dihedral_types_by_element[elements]
		options_by_i = [ [],[],[],[] ]


		if elements in dihedral_types_by_element:
			print elements
			for a in d.atoms:
				a.index2_options.append( set() )
			for t in dihedral_types_by_element[elements]:
				#print '\t', t.index2s
				for i in range(4):
					d.atoms[i].index2_options[-1].add( t.index2s[i] )
		else:
			print 'Error: dihedral', elements, 'does not exist in OPLS file', parameter_file

	for a in atoms:
		print a.element
		for option in a.index2_options:
			print '\t', option

		options = a.index2_options[0]

		for i in xrange(1,len(a.index2_options)):
			options = options.intersection( a.index2_options[i] )

		print '\t\t', options


def opt_opls(molecule, parameter_file=opls_path, taboo_time=100):
	elements, atom_types, bond_types, angle_types, dihedral_types = files.read_opls_parameters(parameter_file)
	atoms, bonds, angles, dihedrals = files.read_cml(molecule, parameter_file=None)

	bond_types_by_index2 = dict( [ (tuple(t.index2s),t) for t in bond_types ] + [ (tuple(reversed(t.index2s)),t) for t in bond_types ] )
	angle_types_by_index2 = dict( [ (tuple(t.index2s),t) for t in angle_types ] + [ (tuple(reversed(t.index2s)),t) for t in angle_types ] )
	dihedral_types_by_index2 = dict( [ (tuple(t.index2s),t) for t in dihedral_types ] + [ (tuple(reversed(t.index2s)),t) for t in dihedral_types ] )

	charges_by_index = dict( [ (t.index,t.charge) for t in atom_types ] )

	for a in atoms:
		a.possible_types = set()
		for t in atom_types:
			if elem_i2s(t.element) == a.element and t.bond_count==len(a.bonded):
				a.possible_types.add(t.index2)
		a.possible_types = list(a.possible_types)
		#a.possible_types.append(0)

	def count_conflicts(types):
		for i,a in enumerate(atoms):
			a.index2 = types[i]
		conflicts = 0
		for b in bonds:
			index2s = (b.atoms[0].index2, b.atoms[1].index2)
			if not index2s in bond_types_by_index2:
				conflicts += 1

		for a in angles:
			index2s = (a.atoms[0].index2, a.atoms[1].index2, a.atoms[2].index2)
			if not index2s in angle_types_by_index2:
				conflicts += 1

		for d in dihedrals:
			index2s_0 = (d.atoms[0].index2, d.atoms[1].index2, d.atoms[2].index2, d.atoms[3].index2)
			index2s_1 = (0,                 d.atoms[1].index2, d.atoms[2].index2, d.atoms[3].index2)
			index2s_2 = (d.atoms[0].index2, d.atoms[1].index2, d.atoms[2].index2,        0)
			in0 = index2s_0 in dihedral_types_by_index2
			in1 = index2s_1 in dihedral_types_by_index2
			in2 = index2s_2 in dihedral_types_by_index2
			if not in0 and not in1 and not in2:
				conflicts += 1
		return conflicts

	import random
	types = [random.choice(a.possible_types) for a in atoms]
	taboo = [0 for a in atoms]
	best = count_conflicts(types)

	step = 0
	for step in range(100000):
		i = random.randint( 0, len(types)-1 )
		for guess in types:
			if taboo[i]>0:
				i = random.randint( 0, len(types)-1 )
			else:
				break
		old_type = types[i]
		types[i] = random.choice(atoms[i].possible_types)

		conflicts = count_conflicts(types)
		if conflicts <= best:
			best = conflicts
			taboo[i] = taboo_time
		else:
			types[i] = old_type

		taboo = [t-1 if t>0 else 0 for t in taboo]

		if step % 10000 == 0:
			print best, conflicts, types
		step += 1

	def types_from_index2(x):
		return [t for t in atom_types if t.index2==x and t.index<=440]

	for i,tt in enumerate( [ types_from_index2(x) for x in types] ):
		#print i, atoms[i].element
		atoms[i].index_options = [t.index for t in tt]
		#for t in tt:
		#	print '\t', t.index, t.notes

	def net_charge(types):
		charge = 0.0
		for t in types:
			charge += charges_by_index[t]
		return charge

	types = [random.choice(a.index_options) for a in atoms]
	taboo = [0 for a in atoms]
	best = net_charge(types)

	for step in range(100000):
		i = random.randint( 0, len(types)-1 )
		for guess in types:
			if taboo[i]>0:
				i = random.randint( 0, len(types)-1 )
			else:
				break
		old_type = types[i]
		types[i] = random.choice(atoms[i].index_options)

		charge = net_charge(types)
		if abs(charge) <= abs(best):
			best = charge
			taboo[i] = taboo_time
		else:
			types[i] = old_type

		taboo = [t-1 if t>0 else 0 for t in taboo]

		if step % 10000 == 0:
			print best, charge, types
		step += 1

	for t in types:
		for t2 in atom_types:
			if t2.index==t:
				print t2.element, t2.notes

def spaced_print(sOut,delim=['\t',' '],buf=4):
	s_len = []
	if type(sOut) == str: sOut = sOut.split('\n')
	if type(delim) == list: delim = ''.join([d+'|' for d in delim])[:-1]
	# Get the longest length in the column
	for i,s in enumerate(sOut):
		s = re.split(delim,s)
		for j,sss in enumerate(s):
			ss = strip_colour(sss)
			try: s_len[j] = len(ss) if len(ss)>s_len[j] else s_len[j] # This makes the part of the list for each column the longest length
			except: s_len.append(len(ss)) # If we are creating a new column this happens
	for i in range(len(s_len)): s_len[i] += buf # Now we add a buffer to each column

	# Compile string output
	for i,s in enumerate(sOut):
		s = re.split(delim,s)
		for j,ss in enumerate(s): s[j] = ss + ''.join([' ']*(s_len[j]-len(strip_colour(ss))))
		sOut[i] = ''.join(s)

	return '\n'.join(sOut)


# http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector/25709323#25709323
def rotation_matrix(axis, theta, units="deg"):
	from numpy import cross, eye, radians
	from scipy.linalg import expm, norm
	if "deg" in units.lower():
		theta = radians(theta)
	return expm(cross(eye(3), axis/norm(axis)*theta))


# Pass a list of atoms, return a list of lists of atoms
def rotate_xyz(frame, theta_0=0, theta_n=360, dt=1, axis=[0,0,1], com=True, last=False):
	from numpy import dot, array, arange, radians
	from scipy.linalg import block_diag
	frames,image = [], array([[a.x,a.y,a.z] for a in frame]).flatten()
	elements = [a.element for a in frame]
	natoms = len(frame)

	if com:
		image = image.reshape((-1,3))
		translate = image.sum(axis=0) / float(natoms)
		image -= translate
		image = image.flatten()

	theta = theta_n if last else theta_0
	while theta <= theta_n:
		r = rotation_matrix(axis, radians(theta))
		R = block_diag(*[r for i in range(natoms)])
		rotated = dot(R,image).reshape((-1,3))
		if com:
			rotated += translate
		frames.append([Atom(e,a[0],a[1],a[2]) for e,a in zip(elements, rotated)])
		theta += dt

	if last: return frames[0]
	return frames


# http://stackoverflow.com/questions/14016898/port-matlab-bounding-ellipsoid-code-to-python/14025140#14025140
# TODO: Proof read later
def mvee(points, tol = 0.001):
    """
    Find the minimum volume ellipse.
    Return A, c where the equation for the ellipse given in "center form" is
    (x-c).T * A * (x-c) = 1
    """
    from numpy import asmatrix, column_stack, argmax, diag, asarray, squeeze, ones
    from numpy.linalg import norm, inv
    elements = [a.element for a in points]
    points = Molecule(points).flatten().reshape((-1,3))
    points = np.asmatrix(points)
    N, d = points.shape
    Q = column_stack((points, ones(N))).T
    err = tol+1.0
    u = ones(N)/N
    while err > tol:
        # assert u.sum() == 1 # invariant
        X = Q * diag(u) * Q.T
        M = diag(Q.T * inv(X) * Q)
        jdx = argmax(M)
        step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
        new_u = (1-step_size)*u
        new_u[jdx] += step_size
        err = norm(new_u-u)
        u = new_u
    c = u*points
    A = inv(points.T*diag(u)*points - c.T*c)/d
    points = asarray(A).flatten().reshape((-1,3))
    centroid = squeeze(asarray(c))
    return points, centroid

# Given a list of atom objects, this will (1) use mvee to generate a minimum centroid around
# the structure, and (2) rotate the ellipsoid and positions to align along the x-axis
def align_centroid(points):
	from scipy.linalg import block_diag
	# Get points and the ellipsoid
	A, centroid = mvee(points)
	points = Molecule(points).flatten().reshape((-1,3))
	npoints = float(len(points))

	# Rotate the ellipsoid
	omega, R = np.linalg.eigh(A)
	A = np.dot(A,R)

	# Rotate the points
	rotation = block_diag(*[R for a in points])
	points = points.flatten()
	points = np.dot(points, rotation).reshape((-1,3))

	# Recenter the points
	molec = Molecule([])
	molec.set_positions(points, new_atom_list=True)
	com = np.array(molec.get_com()) * -1.0
	molec.translate(com)

	return molec.atoms, A


def pysub(job_name, nprocs="1", queue="batch", xhost=None, path=os.getcwd(), remove_nbs=False):
	if type(xhost) is str:
		xhost = [xhost]
	if ".py" in job_name: job_name = job_name.split(".py")[0]
	xhosts = ""
	if xhost is not None:
		xhosts = "##NBS-xhost: " + ", ".join( map(lambda x: '"' + x + '"', xhost) )

	# Setup nbs script
	NBS = '''##NBS-name: "$JOB_NAME$"
##NBS-nproc: $NPROCS$
##NBS-queue: "$QUEUE$"
$XHOST$
source /fs/home/$USER/.zshrc

/fs/home/$USER/anaconda/bin/python2.7 -u $PY_NAME1$.py >> $PY_NAME2$.log 2>&1
'''

	NBS = NBS.replace("$JOB_NAME$",job_name)
	NBS = NBS.replace("$NPROCS$",str(nprocs))
	NBS = NBS.replace("$QUEUE$",queue)
	NBS = NBS.replace("$PY_NAME1$",path + '/' + job_name)
	NBS = NBS.replace("$PY_NAME2$",path + '/' + job_name)
	NBS = NBS.replace("$XHOST$",xhosts)

	NBS_fptr = open(job_name+".nbs",'w')
	NBS_fptr.write(NBS)
	NBS_fptr.close()

	# Submit job
	os.system('jsub ' + job_name + '.nbs')

	if remove_nbs:
		os.system('rm ' + job_name + '.nbs')

def get_pdf(frames, start=0.0, stop=5.0, step=0.1, cutoff=10.0, rho=1.0, quanta=0.001, output=None, persist=False):
	# If passed frames and not an xyz file name, write to xyz
	append = str(int(random.random()*1E12))
	if type(frames) is not str:
		files.write_xyz(frames, "tmp_for_pdf_%s" % append)
		file_name = "tmp_for_pdf_%s" % append
	else:
		file_name = frames

	# Else, we want to ensure file_name is correct
	if file_name.endswith(".xyz"):
		file_name = file_name.split(".xyz")[0]
	if output is None:
		output = file_name

	if stop > cutoff:
		raise Exception("Stopping position should be larger less than or equal to the cutoff.")

	# Make command for debyer
	cmd = "debyer --cutoff=%.2f --quanta=%.2f -g -f%.2f -t%.2f -s%.2f --ro=%.2f -o %s.g %s.xyz" % (cutoff, quanta, start, stop, step, rho, output, file_name)

	# Run debyer and read in the pdf
	os.system(cmd)
	fptr_pdf = open("%s.g" % output, 'r').read().split("\n")
	i = 0
	while fptr_pdf[i].strip().startswith("#"): i += 1
	j = len(fptr_pdf)-1
	while fptr_pdf[j].strip() == "": j -= 1
	fptr_pdf = fptr_pdf[i:j+1]

	pdf = [(float(a.split()[0]), float(a.split()[1])) for a in fptr_pdf]

	if not persist:
		os.system("rm %s.g" % output)
		os.system("rm %s.xyz" % file_name)

	return pdf

def clean_up_folder(path, files_to_remove=[], remove_empty_folders=False, verbose=False):
	if len(files_to_remove) == 0 and remove_empty_folders is False:
		raise Exception("clean_up_folders requires either a file identifier to remove files or the remove_empty_folders flag to be set to True.")
	if path[0] != "/":
		raise Exception("For safety reasons, we require a full path to be used in clean_up_folders.")

	if path.endswith("/"): path = path[:-1]

	if verbose: print("\n---------------------\nCleaning up folder %s" % path)
	if len(files_to_remove) > 0:
		# Remove all empty folders
		ids = " -o -name ".join( map(lambda x: '"' + x + '" -delete ', files_to_remove) )
		ids = "-name " + ids
		cmd_delete_files = "find %s -type f %s" % (path, ids)
		if verbose: print("Running command: %s" % cmd_delete_files)
		os.system(cmd_delete_files)

	if remove_empty_folders:
		cmd_delete_folders = "find %s -type d -empty -delete" % path
		if verbose: print("Running command: %s" % cmd_delete_folders)
		os.system(cmd_delete_folders)
	if verbose: print("Done cleaning folder %s\n---------------------\n\n" % path)
