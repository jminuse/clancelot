import numpy as np
from math import fsum

class LJ:
	atoms, gradient, error = None, None, None
	step, RMS = 0, None
	
	def __init__(self, atoms, sigma=1, epsilon=1, rms_coeff=0, e_coeff=1):
		LJ.atoms = atoms
		LJ.step = 0
		LJ.coords_start = []
		LJ.epsilon = np.float128(epsilon)
		LJ.sigma = np.float128(sigma)
		LJ.energy = np.float128('inf')
		LJ.error = np.float128('inf')
		LJ.RMS_coeff = np.float128(rms_coeff)
		LJ.E_coeff = np.float128(e_coeff)
		for a in LJ.atoms:
			LJ.coords_start += [np.float128(a.x), np.float128(a.y), np.float128(a.z)]

	@staticmethod
	def calculate(coords):
		# Update atom coordinates from line array
		coord_count = 0
		for a in LJ.atoms:
			a.x,a.y,a.z = coords[coord_count], coords[coord_count+1], coords[coord_count+2]
			coord_count += 3

		# Set gradients to 0
		LJ.gradient = [np.array([0.,0.,0.]) for i in range(len(LJ.atoms))]

		# Loop through the atoms
		LJ.energy = np.float128(0.)
		energy_hold = []
		for i,aa in enumerate(LJ.atoms):
			a = np.array([aa.x,aa.y,aa.z])
			gx, gy, gz = [], [], []
			for j,bb in enumerate(LJ.atoms):
				# If comparing against same atom, skip
				if i == j: continue
				
				b = np.array([bb.x,bb.y,bb.z])
				dist = np.linalg.norm(a-b)
				dir = (a-b)/dist

				# Equations from http://www.phys.ubbcluj.ro/~tbeu/MD/C2_for.pdf
				calc_F = np.float128(dir * 48.0 * LJ.epsilon / np.power(dist,2) * ( np.power((LJ.sigma/dist),12) - 0.5*np.power((LJ.sigma/dist),6)))
				calc_E = np.float128(4.0 * LJ.epsilon * ( np.power((LJ.sigma/dist),12) - np.power((LJ.sigma/dist),6)))

				gx.append(np.float128(-calc_F[0]))
				gy.append(np.float128(-calc_F[1]))
				gz.append(np.float128(-calc_F[2]))
				energy_hold.append(np.float128(calc_E))

			x, y, z = fsum(gx), fsum(gy), fsum(gz)
			LJ.gradient[i] = np.array([x,y,z])
		LJ.energy = fsum(energy_hold)



		LJ.energy /= np.float128(2.0) # Remove double counts

		RMS = np.sqrt(np.sum(np.array([np.linalg.norm(grad)**2 for grad in LJ.gradient])))
		LJ.error = LJ.RMS_coeff*RMS + LJ.E_coeff*LJ.energy
		LJ.step += 1
		gradient = []
		for pt in LJ.gradient:
			gradient.append(pt[0])
			gradient.append(pt[1])
			gradient.append(pt[2])
		LJ.gradient = np.array(gradient)

	@staticmethod
	def get_error(coords):
		if LJ.error is None:
			LJ.calculate(coords)
		error = LJ.error
		LJ.error = None
		return error

	@staticmethod
	def get_gradient(coords):
		if LJ.gradient is None:
			LJ.calculate(coords)
		gradient = LJ.gradient
		LJ.gradient = None
		return np.array(gradient)