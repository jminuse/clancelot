import sys,  os, re
import g09, orca, files
from merlin import units

dft, u1, u2, scale, step, out_name, comp, neb_force = 'g09', 'Ha', 'Ha', 1.0, 1, 'out', None, None
title, x_label, y_label, x_range, y_range, x_vals = 'Energy Landscape', 'X-Axis', 'Y-Axis', None, None, None

dft_list = [dft,'orca']

help = '''
scanDFT
---------
A command to view the energy landscape over several configurations.
Note, the START STOP range is inclusive on either end.

scanDFT [Sim_Name%%d] START STOP [Options]

    Flag          Default     Description
-help, -h     :            :  Print this help menu
-dft          :    g09     :  Specify what type of dft simulation you want to
                              get the energy landscape of. Other options include
                              'orca'.
-units, -u    :  kcal/mol  :  Specify the units you want the output to be in.
-scale        :    1.0     :  Scale all energies by this value. Applied AFTER
                              unit conversion from simulation units ('Ha') to
                              -units.
-out, -o      :    out     :  Make an output file with this name holding all xyz
                              coordinates of what you're scanning over.
-step         :    1.0     :  Steps to take between your start and stop range (integer)
-c            :            :  Compile multiple energy landscapes on the same graph.
                              Takes three arguments, separated by commas with no spaces:
                                   char,start,stop
                              The character is a unique identifier in the Sim_Name that will
                              be replaced with values from start to stop (inclusive)
-neb          :            :  In NEB calculations, each iteration after the first does not
                              include the first and last energies. Giving this flag and a
                              run name for the first in the NEB list will tack on these energies
                              to the rest of the simulations
    
-title, -t    :            :  Title for the output graph
-lx           :            :  Label for the x-axis
-ly           :            :  Label for the y-axis
-xrange       :            :  Set the x-axis range
-yrange       :            :  Set the y-axis range
-xvals        :            :  Set a custom label for x-axis (comma separated).

ex: scanDFT water_ 1 10
ex: scanDFT water_%d 1 10
ex: scanDFT water%d_opt 1 10
ex: scanDFT water_^_%d 1 10 -c ^,0,4 -dft orca
ex: scanDFT water_^_%d 1 10 -c ^,2,4 -dft orca -neb water_0_0,water_0_10
ex: scanDFT water_opt_%d 1 10 -t "Water Optimization" -xrange 0,5
'''

if '-h' in sys.argv or '-help' in sys.argv or len(sys.argv) < 3:
	print help
	sys.exit()

# READ IN DATA
run_name = sys.argv[1]
start = int(sys.argv[2])
stop = int(sys.argv[3])

if '-dft' in sys.argv:
	dft = sys.argv[sys.argv.index('-dft') + 1].lower()
	if dft not in dft_list:
		print("Error - %s not recognized for dft." % dft)
		sys.exit()

if [s for s in ['-units','-u'] if s in sys.argv]:
	s = '-u' if '-u' in sys.argv else '-units'
	u2 = sys.argv[sys.argv.index(s) + 1]
	if u2 not in constants.ENERGY:
		print("Error - Energy unit not available. Consider using -scale.")
		sys.exit()

if '-scale' in sys.argv:
	scale = float(sys.argv[sys.argv.index('-scale') + 1])

if '-step' in sys.argv:
	step = int(sys.argv[sys.argv.index('-step') + 1])

if [s for s in ['-o','-out'] if s in sys.argv]:
	s = '-o' if '-o' in sys.argv else '-out'
	out_name = sys.argv[sys.argv.index(s) + 1].replace(' ','_')
if len(out_name) < 5 or out_name[-4:] != '.xyz': out_name += '.xyz'

if '-c' in sys.argv:
	comp = sys.argv[sys.argv.index('-c') + 1].split(',')

if '-neb' in sys.argv:
	neb_force = sys.argv[sys.argv.index('-neb') + 1].split(',')

if [s for s in ['-t','-title'] if s in sys.argv]:
	s = '-t' if '-t' in sys.argv else '-title'
	title = sys.argv[sys.argv.index(s) + 1]

if '-lx' in sys.argv:
	x_label = sys.argv[sys.argv.index('-lx') + 1]
if '-ly' in sys.argv:
	y_label = sys.argv[sys.argv.index('-ly') + 1]

if '-xrange' in sys.argv:
	x_range = sys.argv[sys.argv.index('-xrange') + 1].split(',')
	x_range = [float(x) for x in x_range]
if '-yrange' in sys.argv:
	y_range = sys.argv[sys.argv.index('-yrange') + 1].split(',')
	y_range = [float(y) for y in y_range]

if '-xvals' in sys.argv:
	x_vals = sys.argv[sys.argv.index('-xvals') + 1].split(',')
	x_vals = [float(x) for x in x_vals]


# BEGIN MAKING ENERGY LANDSCAPE
if dft == 'g09':
	read = g09.read
elif dft == 'orca':
	read = orca.read
else:
	print("Error - Cannot proceed with DFT as %s." % dft)
	sys.exit()

first_E, last_E = None, None
first_frame, last_frame = None, None
if neb_force is not None:
	first_E = read(neb_force[0]).energies[-1]
	first_frame = read(neb_force[0]).atoms
	last_E = read(neb_force[1]).energies[-1]
	last_frame = read(neb_force[1]).atoms

energies, frames = [], []
# Loop through energies
if comp == None: comp = [None, 0, 0]
for c in range(int(comp[1]), int(comp[2])+1):
	run_hold = run_name.replace(comp[0],str(c)) if comp[0] is not None else run_name
	tmp_E, tmp_frames = [], []

	if neb_force is not None:
		tmp_frames.append(first_frame)
		tmp_E.append(first_E)

	for i in range(start, stop+1, step):
		# Get run name for this iteration
		run = run_hold+str(i) if run_hold.find('%')==-1 else run_hold % i
		data = read(run)
		tmp_E.append(data.energies[-1])
		tmp_frames.append(data.atoms)

	if neb_force is not None:
		tmp_frames.append(last_frame)
		tmp_E.append(last_E)

	energies.append(tmp_E)
	frames = tmp_frames

# Adjust energies
E_offset = energies[0][0]
for i in range(len(energies)):
	for j,e in enumerate(energies[i]):
		energies[i][j] = units.convert_energy(u1,u2,e-E_offset) * scale

# Plot energies
def plot(yy,start_val,x_label,y_label,title,x_range,y_range):
	import matplotlib.pyplot as plt

	low_x, high_x = start_val,len(yy[0])-1
	low_y, high_y = float('inf'), float('-inf')
	if x_vals is None:
		for i,y in enumerate(yy):
			plt.plot(y,marker='.',label=str(int(start_val) + i))
			if min(y) < low_y: low_y = min(y)
			if max(y) > high_y: high_y = max(y)
	else:
		low_x = min(x_vals)
		high_x = max(x_vals)
		for y in yy:
			plt.plot(x_vals,y,marker='.',label=str(int(start_val) + i))
			if min(y) < low_y: low_y = min(y)
			if max(y) > high_y: high_y = max(y)

	plt.xlabel(x_label)
	plt.ylabel('%s (%s)' % (y_label,u2))
	plt.title(title)

	if x_range is None: x_range = [low_x, high_x]
	if y_range is None: y_range = [low_y, high_y*1.05]

	plt.axis([x_range[0], x_range[1], y_range[0], y_range[1]])

	plt.legend()
	plt.show()

plot(energies,start,x_label,y_label,title,x_range,y_range)

# Write files
files.write_xyz(frames,out_name[:-4])
