import os, string, sys, re, shutil, copy
from subprocess import Popen
import utils, units, log, files

def job(run_name, route, atoms=[], extra_section='', queue='batch', procs=1, charge_and_multiplicity='0 1', title='run by gaussian.py', blurb=None, eRec=True, force=False, previous=None, neb=[False,None,None,None], err=False, lint=False):
	os.chdir('orca')
	f = open(run_name+'.orca', 'w')
	f.write(route+'\n*xyz '+charge_and_multiplicity+'\n')
	for a in atoms:
		f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
	f.write('*\n')
	os.system('/fs/home/jms875/build/orca/orca_3_0_2_linux_x86-64/orca %s.orca > %s.out &' % (run_name,run_name))
	os.chdir('..')
