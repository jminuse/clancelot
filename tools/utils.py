import os, sys, math, copy, subprocess, time
import numpy
import files

class Struct:
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)
	def __repr__(self):
		return str( dict([ (a,None) if type(self.__dict__[a]) in (list,dict) else (a,self.__dict__[a]) for a in self.__dict__]) )

class Atom():
	def __init__(self, element, x, y, z, index=None, type=None, molecule_index=1):
		self.element = element
		self.x = x
		self.y = y
		self.z = z
		self.index = index
		self.type = type
		self.molecule_index = molecule_index

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

def get_bonds(atoms):
	bonds = []
	for i,a in enumerate(atoms):
		a.bonded = []
		for b in atoms[i+1:]:
			dd = dist_squared(a,b)
			if (a.element not in [1,'H'] and b.element not in [1,'H'] and dd<2**2) or (dd < 1.2**2 and (a.element in [1,'H'])!=(b.element in [1,'H']) ) or (dd < 2.8**2 and (a.element in ['Pb',82] or b.element in ['Pb',82]) ):
				bonds.append( Bond(a,b) ) #offset from current, distance
	return bonds

def get_angles_and_dihedrals(atoms):
	angles = [];
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
	while low < high:
		yield low
		low += step

def quat_to_mat(q): #quat = [w i j k]
	d,b,c,a = q
	return [ [a**2+b**2-c**2-d**2, 2*b*c-2*a*d, 2*b*d+2*a*c],
	[2*b*c+2*a*d, a**2-b**2+c**2-d**2, 2*c*d-2*a*b],
	[2*b*d-2*a*c, 2*c*d-2*a*b, a**2-b**2-c**2+d**2]
	]

def matvec(m,v):
	return (m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2], m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2], m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2])
	
def matmat(a,b):
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

elements_by_atomic_number = ['','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Fl','Uup','Lv','Uus','Uuo']


class Molecule():	
	def __init__(self, atoms_or_filename_or_all, bonds=None, angles=None, dihedrals=None, parameter_file='oplsaa.prm', extra_parameters={}): #set atoms, bonds, etc, or assume 'atoms' contains all those things if only one parameter is passed in
		if type(atoms_or_filename_or_all)==type('string'):
			self.filename = atoms_or_filename_or_all
			atoms, bonds, angles, dihedrals = files.read_cml(self.filename, parameter_file=parameter_file, extra_parameters=extra_parameters)
		elif not bonds:
			atoms, bonds, angles, dihedrals = atoms_or_filename_or_all
		else:
			atoms = atoms_or_filename_or_all
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


class System():
	def __init__(self, name=None, box_size=None):
		self.atoms, self.bonds, self.angles, self.dihedrals = [], [], [], []
		self.box_size = box_size
		self.name = name
		self.molecules = []

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
		new_molecule.bonds = self.atoms[-len(molecule.bonds):]
		new_molecule.angles = self.atoms[-len(molecule.angles):]
		new_molecule.dihedrals = self.atoms[-len(molecule.dihedrals):]
		self.molecules.append( new_molecule )

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

def procrustes(frames, count_atoms=None):
	if not count_atoms: count_atoms = range(len(frames[0]))
	for s in frames:
		center_x = sum([a.x for i,a in enumerate(s) if i in count_atoms])/len(count_atoms)
		center_y = sum([a.y for i,a in enumerate(s) if i in count_atoms])/len(count_atoms)
		center_z = sum([a.z for i,a in enumerate(s) if i in count_atoms])/len(count_atoms)
		for a in s:
			a.x -= center_x
			a.y -= center_y
			a.z -= center_z
	#rotate all frames to be as similar to their neighbors as possible
	from scipy.linalg import orthogonal_procrustes
	for i in range(1,len(frames)): #rotate all frames to optimal alignment
		#only count spring-held atoms for finding alignment
		count_atoms_1 = [(a.x,a.y,a.z) for j,a in enumerate(frames[i]) if j in count_atoms]
		count_atoms_2 = [(a.x,a.y,a.z) for j,a in enumerate(frames[i-1]) if j in count_atoms]
		rotation = orthogonal_procrustes(count_atoms_1,count_atoms_2)[0]
		#rotate all atoms into alignment
		for a in frames[i]:
			a.x,a.y,a.z = matvec(rotation, (a.x,a.y,a.z))

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

def pretty_xyz(name,R_MAX=1,F_MIN=1,F_MAX=50,CENTER=None,outName=None,write_xyz=False,verbose=False):
	# Get data as either frames or a file
	if type(name)==type(''): frames = files.read_xyz(name)
	elif type(name)==type([]): frames = name
	else:
		print "Error - Invalid name input.  Should be either the name of an xyz file or a list.", sys.exc_info()[0]
		exit()

	# Loop till we're below R_MAX
	while 1:
		# Check if we're done
		r2 = max(motion_per_frame(frames))
		if r2 < R_MAX: break

		if len(frames) > F_MAX:
			print "-------------------------------------------------------"
			print motion_per_frame(frames)
			print "-------------------------------------------------------"
			print "\n\nError - Could not lower motion below %lg in %d frames." % (R_MAX,F_MAX), sys.exc_info()[0]
			exit()
		else:
			if verbose: print "Currently Frames = %d\tr2 = %lg" % (len(frames),r2)

		# If not, find largest motion_per_frame
		index, maxVal = 0, 0
		tmp = motion_per_frame(frames)
		for i,t in enumerate(tmp):
			if t > maxVal:
				index, maxVal = i, t
		i = index
		# Now, split the list, interpolate, and regenerate
		if i>0 and i < len(frames) - 3:
			f_low = frames[:i-1]
			f_high = frames[i+2:]
			f_mid = interpolate(frames[i-1],frames[i+2],4)
			frames = f_low + f_mid + f_high
		elif i == 0:
			f_high = frames[i+2:]
			f_mid = interpolate(frames[i],frames[i+1],4)
			frames = f_mid + f_high
		else:
			f_low = frames[:i-1]
			f_mid = interpolate(frames[i-1],frames[i],4)
			frames = f_low + f_mid

		if verbose: print "\tInterpolated %d,%d ... %lg" % (index-1,index+1,max(motion_per_frame(frames)))

	if CENTER != None: center_frames(frames,CENTER)

	if write_xyz: files.write_xyz(frames,'pretty_xyz' if outName==None else outName)
	else: return frames