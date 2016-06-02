import os, sys
import math, copy, subprocess, time, numpy, re
import files, constants
from units import elem_i2s
from warning import warn

class Struct:
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
	def __repr__(self):
		return str( dict([ (a,None) if type(self.__dict__[a]) in (list,dict) else (a,self.__dict__[a]) for a in self.__dict__]) )

class Atom():
	def __init__(self, element, x, y, z, index=None, type=None, molecule_index=1, bonded=[], type_index=None):
		self.element = element
		self.x = x
		self.y = y
		self.z = z
		self.index = index
		self.molecule_index = molecule_index
		self.bonded=bonded
		self.type_index=type_index
		# When parameterized by OPLS, 'type' dict contains: {'bond_count': 3, 'index': 588, 'notes': '1,10-Phenanthroline C2', 'vdw_r': 3.55, 'element': 6, 'vdw_e': 0.07, 'charge': 0.392, 'mass': 12.011, 'index2': 48, 'element_name': 'CA'}
		# Object May also contain lammps_type, mass, and charge
		self.type = type

	def translate(self, v):
		self.x+=v[0]; self.y+=v[1]; self.z+=v[2]

	def printSelf(self):
		if self.type_index:
			text = '%s, (%3.3f, %3.3f, %3.3f), index: %d, type_index: %d\n' % (self.element, self.x, self.y, self.z, self.index, self.type_index)
		else:
			text = '%s, (%3.3f, %3.3f, %3.3f), %d\n' % (self.element, self.x, self.y, self.z, self.index)

		if self.type:
			text += str(self.type)
			text += '\n'
		return text

	def __str__(self):
		return self.printSelf()

class Bond():
	def __init__(self, a, b, type=None):
		self.atoms = (a,b)
		self.type = type

class Angle():
	def __init__(self, a, b, c, type=None):
		self.atoms = (a,b,c)
		self.type = type

class Dihedral():
	def __init__(self, a, b, c, d, type=None):
		self.atoms = (a,b,c,d)
		self.type = type

class Job(): #a job on the queue
	def __init__(self, name):
		self.name = name
	def wait(self):
		while True:
			jlist = subprocess.Popen('jlist', shell=True, stdout=subprocess.PIPE).communicate()[0]
			if (' '+self.name+' ') in jlist:
				time.sleep(60)
			else: break

# A generic class to hold dft data
class DFT_out():
	def __init__(self, name, dft='g09'):
		self.name = name
		self.dft = dft.lower()

		# Initialize everything as None
		self.route = None
		self.frames = None
		self.atoms = None
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
class sim_out():
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
			if (a.element not in [1,'H'] and b.element not in [1,'H'] and dd<2**2) or (dd < 1.2**2 and (a.element in [1,'H'])!=(b.element in [1,'H']) ) or (dd < 2.8**2 and (a.element in ['Pb',82] or b.element in ['Pb',82]) ):
				bonds.append( Bond(a,b) )
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
				angles.append( Angle(a,center,b) )
	dihedral_set = {}
	for angle in angles:
		for a in angle.atoms[0].bonded:
			if a is angle.atoms[1]: continue
			dihedral = (a,) + angle.atoms
			if tuple(reversed(dihedral)) not in dihedral_set:
				dihedral_set[dihedral] = True

		for b in angle.atoms[2].bonded:
			if b is angle.atoms[1]: continue
			dihedral = angle.atoms + (b,)
			if tuple(reversed(dihedral)) not in dihedral_set:
				dihedral_set[dihedral] = True
	dihedrals = [Dihedral(*d) for d in dihedral_set.keys()]

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
	import random
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

class Molecule():
	def __init__(self, atoms_or_filename_or_all, bonds=None, angles=None, dihedrals=None, parameter_file='oplsaa.prm', extra_parameters={}, test_charges=True, allow_errors=False): #set atoms, bonds, etc, or assume 'atoms' contains all those things if only one parameter is passed in
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
			a.x, a.y, a.z = matvec(m, (a.x, a.y, a.z))
	def translate(self, v):
		for a in self.atoms:
			a.x+=v[0]; a.y+=v[1]; a.z += v[2]
	def rand_rotate(self):
		rand_m = rand_rotation()
		self.rotate(rand_m)

	# Print all atoms
	def printAtoms(self):
		text = ''
		for atom in self.atoms:
			text += atom.printSelf()
		return text

	# When printing molecule, print all atoms
	def __str__(self):
		return self.printAtoms()


class System():
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

	def add(self, molecule, x=0.0, y=0.0, z=0.0):
		atom_offset = len(self.atoms)
		for a in molecule.atoms:
			new_atom = copy.copy(a)
			new_atom.index=a.index+atom_offset
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

	# Print all atoms
	def printAtoms(self):
		text = ''
		for atom in self.atoms:
			text += atom.printSelf()
		return text

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
	cache_key = a.x+b.x+c.x+d.x

	sqrt = math.sqrt
	cos = math.cos

	#a,b,c,d = d,c,b,a

	vb1x, vb1y, vb1z = a.x-b.x, a.y-b.y, a.z-b.z
	vb2x, vb2y, vb2z = c.x-b.x, c.y-b.y, c.z-b.z
	vb3x, vb3y, vb3z = d.x-c.x, d.y-c.y, d.z-c.z

	# c0 calculation

	sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
	sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
	sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

	rb1 = sqrt(sb1)
	rb3 = sqrt(sb3)

	c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3

	# 1st and 2nd angle

	b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z
	b1mag = sqrt(b1mag2)
	b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z
	b2mag = sqrt(b2mag2)
	b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z
	b3mag = sqrt(b3mag2)

	ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z
	r12c1 = 1.0 / (b1mag*b2mag)
	c1mag = ctmp * r12c1

	ctmp = vb2x*vb3x + vb2y*vb3y + vb2z*vb3z
	r12c2 = 1.0 / (b2mag*b3mag)
	c2mag = -ctmp * r12c2

	# cos and sin of 2 angles and final c

	sin2 = max(1.0 - c1mag*c1mag,0.0)
	sc1 = sqrt(sin2)
	#if sc1 < SMALL: sc1 = SMALL
	sc1 = 1.0/sc1

	sin2 = max(1.0 - c2mag*c2mag,0.0)
	sc2 = sqrt(sin2)
	#if sc2 < SMALL: sc2 = SMALL
	sc2 = 1.0/sc2

	s1 = sc1 * sc1
	s2 = sc2 * sc2
	s12 = sc1 * sc2
	c = (c0 + c1mag*c2mag) * s12

	cx = vb1y*vb2z - vb1z*vb2y
	cy = vb1z*vb2x - vb1x*vb2z
	cz = vb1x*vb2y - vb1y*vb2x
	cmag = sqrt(cx*cx + cy*cy + cz*cz)
	dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag


	if c>1.0: c = 1.0
	if c<-1.0: c = -1.0

	phi = math.acos(c)
	if dx < 0.0:
		phi *= -1.0
	phi *= -1.0

	return phi, math.cos(phi), math.cos(2*phi), math.cos(3*phi), math.cos(4*phi)

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

def pretty_xyz(name,R_MAX=1,F_MAX=50,PROCRUSTS=False,outName=None,write_xyz=False,verbose=False):
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
		if PROCRUSTS: procrustes(frames)
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

	if PROCRUSTS: procrustes(frames)

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

def opls_options(molecule, parameter_file='oplsaa.prm'):
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


def opt_opls(molecule, parameter_file='oplsaa.prm', taboo_time=100):
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
		for j,ss in enumerate(s):
			try: s_len[j] = len(ss) if len(ss)>s_len[j] else s_len[j] # This makes the part of the list for each column the longest length
			except: s_len.append(len(ss)) # If we are creating a new column this happens
	for i in range(len(s_len)): s_len[i] += buf # Now we add a buffer to each column

	# Compile string output
	for i,s in enumerate(sOut):
		s = re.split(delim,s)
		for j,ss in enumerate(s): s[j] = ss + ''.join([' ']*(s_len[j]-len(ss)))
		sOut[i] = ''.join(s)

	return '\n'.join(sOut)

