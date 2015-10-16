from merlin import *

def job(run_name, route, atoms=[], extra_section='', queue=None, procs=1, charge_and_multiplicity='0 1', title='', blurb=None, force=False, previous=None, neb=[False,None,None,None], lint=False):
	# Generate the orca input file
	os.chdir('orca')

	# If running on system with more than one core, tell orca
	if procs > 1:
		route += ' PAL%d' % procs

	# Get input for orca formatted correctly
	inp = route.strip()+'\n'+extra_section.strip()+'\n*xyz '+charge_and_multiplicity+'\n'
	for a in atoms:
		inp += '%s %f %f %f\n' % (a.element, a.x, a.y, a.z)
	inp += '*\n'

	# Write the orca file
	f = open(run_name+'.orca', 'w')
	f.write(inp)
	f.close()

	# Run the simulation
	if queue is None:
		os.system('/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca %s.orca > %s.out &' % (run_name,run_name))
	else:
		NBS = '''#!/bin/bash
##NBS-name: "%s"
##NBS-nproc: %d
##NBS-queue: "%s"

export PATH=/fs/europa/g_pc/ompi_1_6_5/bin:/fs/europa/g_pc/orca_3_0_3_linux_x86-64:$PATH
export LD_LIBRARY_PATH=/fs/europa/g_pc/ompi_1_6_5/lib:$LD_LIBRARY_PATH

/fs/europa/g_pc/orca_3_0_3_linux_x86-64/orca %s.orca >> %s.out 2>&1
''' % (run_name, procs, queue, run_name, run_name)
		f = open(run_name+'.nbs', 'w')
		f.write(NBS)
		f.close()
		os.system('jsub %s.nbs' % run_name)

	# Return to the appropriate directory
	os.chdir('..')