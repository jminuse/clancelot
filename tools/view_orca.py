import os, sys
import files, orca

# Get input
if len(sys.argv) < 1:
	print("Error - viewo requires input.")
	sys.exit()
input = sys.argv[1]

# Check we we want to view output in vmd (default = True)
vmd = True
if '-vmd' in sys.argv:
	tmp = sys.argv[sys.argv.index('-vmd')+1]
	if tmp.lower() in ['false','no','n','f','0']:
		vmd = False

# Check we we want to view energies (default = True)
ener = True
if '-ener' in sys.argv:
	tmp = sys.argv[sys.argv.index('-ener')+1]
	if tmp.lower() in ['false','no','n','f','0']:
		ener = False

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
if not os.path.isfile('orca/%s/%s.orca.trj' % (input,input)):
	print("Error - Trajectory file orca/%s/%s.orca.trj does not exist." % (input,input))
	sys.exit()
os.system('cp orca/%s/%s.orca.trj %s' % (input, input, name))

# Read in data from file
energies = orca.parse_atoms(input,get_atoms=False,parse_all=True)[0]
if type(energies) != list: energies = [energies]

if len(energies) < 1:
	print("No data available.")
	sys.exit()

if ener:
	print('\n'.join(str(e) for e in energies))

if vmd:
	os.system('/fs/europa/g_pc/vmd-1.9 '+name+' > /dev/null')
