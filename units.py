from constants import ENERGY

def convert_energy(e0, e1, e_val):
	if e0 == e1: return e_val
	val = e_val * ENERGY[e0] # This many joules
	return val/ENERGY[e1] # This many of unit e1