import os, sys
import subprocess
import lammps_job, lammps_log
import units, utils, sysconst, constants, files
from getpass import getuser

USERNAME = getuser()

# One can easily change defaults here (if they do, probably change the
# help text below accordingly).
md, u1, u2, scale, out_name, vmd = 'lammps', 'kcal/mol', 'kcal/mol', 1.0, 'final.lammpstrj', False
md_list = [md]

if '-h' in sys.argv or '-help' in sys.argv:
	print('''
chkMD
---------
A command to quickly get a glimpse of a MD simulation.
chkMD [Sim_Name] [Options]

    Flag          Default              Description
-md         	:  lammps  			:  Specify what type of md simulation you want to
                          				parse. Currently only lammps supported
-help, -h  		:             		:  Print this help menu
-units, -u 		:  kcal/mol   		:  Specify the units you want the output to be in.
                                		By default this is kcal/mol for lammps units real.
-scale      	:  1.0        		:  Scale all energies by this value
-out, -o    	:  final.lammpstrj  :  Select a particular dump file to use for output
-vmd, -v    	:                	:  Opens select dump file in vmd. Flag turns on.

ex. chkMD water -u kT_300
''')
	sys.exit()

# Get simulation name
run_name = sys.argv[1]

# Get the md type
if '-md' in sys.argv:
	md = sys.argv[sys.argv.index('-md') + 1].lower()
	if md not in md_list:
		print("Error - %s not recognized for md." % md)
		sys.exit()

# Get units
if '-u' in sys.argv or '-units' in sys.argv:
	s = '-u' if '-u' in sys.argv else '-units'
	u2 = sys.argv[sys.argv.index(s) + 1]
	if u2 not in constants.ENERGY:
		print("Error - Energy unit not available. Consider using -scale.")
		sys.exit()

# Get scale
if '-scale' in sys.argv:
	scale = float(sys.argv[sys.argv.index('-scale') + 1])

# Get dump file name
if '-out' in sys.argv or '-o' in sys.argv:
	s = '-o' if '-o' in sys.argv else '-out'
	out_name = sys.argv[sys.argv.index(s) + 1].replace(' ','_')

# Get VMD display status
if '-vmd' in sys.argv or '-v' in sys.argv:
	vmd = True

# Read in data
if md == 'lammps':
	try:
		lg, data_trj = lammps_job.read(run_name, read_atoms=False, read_timesteps=True, read_num_atoms=False, read_box_bounds=False)
	except IOError:
		print("Error - lammps simulation %s does not exist." % run_name)
		sys.exit()
else:
	print("MD type %s not available..." % md)
	sys.exit()

# Get the header information
head = 'Job Name: %s\n' % run_name
head += 'MD calculation via %s\n' % md
body, tail = '', ''

# Check how many timesteps were written to the dump file and what the last timestep is
try:
	body += '# of Timesteps: %d, Final Timestep: %d\n' % (len(data_trj.timesteps), data_trj.final_timestep)
except TypeError:
	print("No atomic coordinates available yet...")
except:
	print("An unexpected error has occurred.")
	sys.exit()

# Write when the files were last written to
if lg.last_modified != 'Null':
	tail += 'Log file last modified: %s\n' % (lg.last_modified)
if data_trj.last_modified != 'Null':
	tail += 'Dump file last modified: %s\n' % (data_trj.last_modified)

length = max([len(tmp) for tmp in head.split('\n')] + [len(tmp) for tmp in body.split('\n')] + [len(tmp) for tmp in tail.split('\n')])
dash = '\n'+''.join(['-']*length)+'\n'

print(dash+head+dash+body+dash+tail)

try:
	if len(data_trj.frames) > 0:
		if vmd:
			os.system('"'+sysconst.vmd_path + '" ' + me + out_name)
except TypeError:
	print("No atomic coordinates available yet...")
except:
	print("An unexpected error has occurred.")
	sys.exit()
