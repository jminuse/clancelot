from merlin import *
from subprocess import Popen
from copy import deepcopy
from shutil import copyfile
import re

# This module provides a framework for submitting lammps jobs sent through the NBS queue

# A function to run an LAMMPS Simulation. Requires a run name and a string of lammps code (run_name and input_script)
def job(run_name, input_script, system, queue=None, procs=1, email=''):
	if len(run_name) > 31 and queue is not None:
		raise Exception("Job name too long (%d) for NBS. Max character length is 31." % len(run_name))

	# Change to correct directory
	os.system('mkdir -p lammps/%s' % run_name)
	os.chdir('lammps/%s' % run_name)

	# Generate the lammps data file
	files.write_lammps_data(system, pair_coeffs_included=True)

	# Write the lammps input script. Expects lines of lammps code
	f = open(run_name+'.in', 'w')
	f.write(input_script)
	f.close()

	# Write email settings if set, otherwise blank out the line
	if len(email) > 0:
		email_text = '##NBS-email: '+email
	else:
		email_text = '# No email selected'

	# Run the simulation
	if queue is None:
		os.system('/fs/home/yma3/Software/lammps_git/src/lmp_g++ -in %s.in -log %s.log' % (run_name, run_name))
		#process_handle = Popen('/fs/home/yma3/Software/lammps_git/src/lmp_g++ -in %s.in -log %s.log > %s.out' % (run_name, run_name, run_name), shell=True)
		print('Finished local lammps job')
	elif queue=='debug':
		print 'Would run', run_name
	else: #adapted from orcasub
		NBS = '''#!/bin/sh
##NBS-unique: yes
##NBS-name: '''+run_name+'''
##NBS-nproc: '''+str(procs)+'''
##NBS-queue: '''+queue+'''
'''+email_text+'''

# Shared LAMMPS Library Configuration
export PYTHONPATH=/fs/home/yma3/Software/lammps_git/python:$PYTHONPATH
export PYTHONPATH=/fs/home/yma3/Projects/Forcefields/OPLS:$PYTHONPATH
export LD_LIBRARY_PATH=/fs/home/yma3/Software/lammps_git/src:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/mpich2/icse/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/fs/home/yma3/usr/local/lib/fftw3/lib:$LD_LIBRARY_PATH

#$NBS_PATH/mpiexec -xlate /usr/common/etc/nbs/mpi.xlate /fs/home/yma3/Software/lammps_git/src/lmp_g++ -partition 1x'''+str(procs)+''' -in ''' + (os.getcwd()+'/'+run_name) + '''.in -log ''' + (os.getcwd()+'/'+run_name) + '''.log > ''' + (os.getcwd()+'/'+run_name) + '''.out

$NBS_PATH/mpiexec -xlate /usr/common/etc/nbs/mpi.xlate /fs/home/yma3/Software/lammps_git/src/lmp_g++ -in ''' + (os.getcwd()+'/'+run_name) + '''.in -echo log -log ''' + (os.getcwd()+'/'+run_name) + '''.log

'''
		f = open(run_name+'.nbs', 'w')
		f.write(NBS)
		f.close()
		os.system('jsub %s.nbs' % run_name)

	# Copy run script
	fname = sys.argv[0]
	if '/' in fname: fname = fname.split('/')[-1]
	try:
		copyfile('../../%s' % fname,fname)
	except IOError:
		# Submitted a job oddly enough that sys.argv[0] is not the original python file name, so don't do this
		pass

	# Return to the appropriate directory
	os.chdir('../..')

	if queue is None:
		return 'finished'
		#return process_handle
	else:
		return utils.Job(run_name)

