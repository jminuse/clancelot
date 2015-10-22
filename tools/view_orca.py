import os, sys
import files, orca

# Get input
if len(sys.argv) < 1:
	print("Error - viewo requires input.")
	sys.exit()
input = sys.argv[1]

# Check we we want to view output in vmd (default = True)
if '-vmd' in sys.argv:
	tmp = sys.argv[sys.argv.index('-vmd')+1]
	if tmp.lower() in ['false','no','n','f','0']:
		vmd = False
	else:
		vmd = True

# Check we we want to view energies (default = True)
if '-ener' in sys.argv:
	tmp = sys.argv[sys.argv.index('-ener')+1]
	if tmp.lower() in ['false','no','n','f','0']:
		ener = False
	else:
		ener = True

# Check if we want to name output (default = out.xyz)
if '-out' in sys.argv:
	name = sys.argv[sys.argv.index('-out')+1]
elif '-o' in sys.argv:
	name = sys.argv[sys.argv.index('-o')+1]
else:
	name = 'out'
if len(name) > 4:
	if name[-4:] != '.xyz':
		name += '.xyz'
else:
	name += '.xyz'

# Get the xyz file
os.system('cp orca/%s/%s.traj %s' % (input, input, name))

# Read in data from file
energies = orca.parse_atoms(input,get_atoms=False,parse_all=True)
if type(energies) != list: energies = [energies]

if len(energies) < 1:
	print("No data available.")
	sys.exit()

if ener:
	print('\n'.join(energies))

if vmd:
	os.system('/fs/europa/g_pc/vmd-1.9 '+name+' > /dev/null')
