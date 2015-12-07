import constants
from constants import ENERGY, DISTANCE, PERIODIC_TABLE

def convert_energy(e0, e1, e_val):
	if e_val == 0: return 0
	if e0 == e1: return e_val
	if len(e0) > 2 and e0[0:2] == 'kT':
		val = e_val * constants.K_b * float(e0[3:])
	else:
		val = e_val * ENERGY[e0] # This many joules
	if len(e1) > 2 and e1[0:2] == 'kT':
		return val/(constants.K_b * float(e1[3:]))
	return val/ENERGY[e1] # This many of unit e1

def convert_dist(d0, d1, d_val):
	if d_val == 0: return 0
	if d0 == d1: return d_val
	val = d_val * DISTANCE[d0] # This many joules
	return val/DISTANCE[d1] # This many of unit e1

def elem_i2s(elem_int):
	try: i = int(elem_int) # Check if it's already a symbol, if so return it
	except: return elem_int
	return PERIODIC_TABLE[i]['sym'] # If not, convert it

def elem_s2i(elem_sym):
	# CHECK IF elem_int IS ALREADY ALPHANUMERIC!

	# Else convert it
	for i in range(1,len(PERIODIC_TABLE)):
		if PERIODIC_TABLE[i]['sym'] == elem_sym:
			return i

	# Return -1 if you couldn't
	return -1

def elem_weight(elem):
	if type(elem) == str: return PERIODIC_TABLE[elem_s2i(elem)]['weight']
	if type(elem) == int: return PERIODIC_TABLE[elem]['weight']
	print("Warning - No weight found for %s!" % str(elem))
	return -1

def convert(old,new,val):
	if val == 0: return 0
	
	a,b = old.split('/')
	a2,b2 = new.split('/')

	if a in ENERGY:	new_val = convert_energy(a,a2,val)
	else: new_val = convert_dist(a,a2,val)

	if new_val != 0: new_val = 1./new_val

	if b in ENERGY: new_val = convert_energy(b,b2,new_val)
	else: new_val = convert_dist(b,b2,new_val)

	return 1./new_val
