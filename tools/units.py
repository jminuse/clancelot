from constants import ENERGY, PERIODIC_TABLE

def convert_energy(e0, e1, e_val):
	if e0 == e1: return e_val
	val = e_val * ENERGY[e0] # This many joules
	return val/ENERGY[e1] # This many of unit e1

def elem_i2s(elem_int):
	try: i = int(elem_int) # Check if it's already a symbol, if so return it
	except: return elem_int
	return PERIODIC_TABLE[i]['sym'] # If not, convert it