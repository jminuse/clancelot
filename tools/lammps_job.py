from merlin import *
from lammps_log import lammps_log
from subprocess import Popen
from copy import deepcopy
from shutil import copyfile
import re
import cPickle as pickle
from warnings import warn
import sys

# This module provides a framework for submitting and reading lammps jobs sent through the NBS queue.
# LAMMPS is inherently flexible. However, in order to provide a smooth user experience, this module expects certain dump files and formats.
# 1) The job will be created in lammps/run_name/
# 2) The log file is named run_name.log
# 3) The data file is named run_name.data
# 4) The default dump file is named dump.run_name.lammpstrj
# 5) The final simulation configuration can be found in final.lammpstrj

# A function to read both a log file containing thermo information and (optional) a lammpstrj file
def read(run_name, trj_file='', read_atoms=True, read_timesteps=True, read_num_atoms=True, read_box_bounds=True):
	# Format log file name
	if run_name.startswith('/'): #allow absolute paths as filenames
		log_path = run_name
	else:
		log_path = 'lammps/%s/%s.log' % (run_name,run_name)

	# Check if log file exists, and open
	if not os.path.isfile(log_path):
		raise IOError('Expected lammps log file does not exist at %s' % (log_path))
		sys.exit()
	else:
		lg = lammps_log(log_path)
		lg.last_modified = files.last_modified(log_path)

	# If no trj_file selected, try default name of dump.run_name.lammpstrj
	if trj_file is not None:
		if trj_file == '':
			trj_path = 'lammps/%s/dump.%s.lammpstrj' % (run_name, run_name)
		# Allow absolute paths as filenames
		elif run_name.startswith('/'):
			trj_path = trj_file
		# Open the specified file in the run_name folder
		else:
			trj_path = 'lammps/%s/%s' % (run_name, trj_file)

		# Try to import lammpstrj file exists
		data_trj = files.read_lammpstrj(trj_path, read_atoms=read_atoms, read_timesteps=read_timesteps, read_num_atoms=read_num_atoms, read_box_bounds=read_box_bounds)
		data_trj.last_modified = files.last_modified(trj_path)

		return lg, data_trj
	else:
		return lg, None

# A function to run an LAMMPS Simulation. Requires a run name and a string of lammps code (run_name and input_script)
def job(run_name, input_script, system, queue=None, procs=1, email='',
		pair_coeffs_included=True, hybrid_pair=False, hybrid_angle=True, TIP4P=False):
	if len(run_name) > 31 and queue is not None:
		raise Exception("Job name too long (%d) for NBS. Max character length is 31." % len(run_name))

	# Change to correct directory
	os.system('mkdir -p lammps/%s' % run_name)
	os.chdir('lammps/%s' % run_name)

	# Generate the lammps data file
	files.write_lammps_data(system, pair_coeffs_included=pair_coeffs_included,
							hybrid_pair=hybrid_pair, hybrid_angle=hybrid_angle)

	# If TIP4P water is used, read the data file for the correct TIP4P types then inject into the input script
	if TIP4P:
		otype, htype, btype, atype = read_TIP4P_types(run_name+'.data')

		# Looking for TIP4P related strings:
		# pair_style lj/cut/tip4p/long -1 -1 -1 -1 0.1250 12.0
		new_text = 'lj/cut/tip4p/long %d %d %d %d' % (otype, htype, btype, atype)
		input_script = input_script.replace('lj/cut/tip4p/long -1 -1 -1 -1', new_text)

		# fix RIGID_WATER all shake 0.0001 20 1000 b -1 a -1
		new_text = 'b %d a %d' % (btype, atype)
		input_script = input_script.replace('b -1 a -1', new_text)

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
		os.system('lammps -in %s.in -log %s.log' % (run_name, run_name))
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

def read_TIP4P_types(data_file):
	otype, htype, btype, atype = -1, -1, -1, -1

	# Iterate file line by line. This reduces the memory required since it only loads the current line
	with open(data_file) as f:
		for line in f:
			a = line.split()

			if len(a) == 3:
				# Check for oxygen
				if a[1] == '0.162750' and a[2] == '3.164350':
					otype = int(a[0])

				# Check for hydrogen
				if a[1] == '0.000000' and a[2] == '0.000000':
					htype = int(a[0])

				# Check for TIP4P bond
				if a[1] == '7777.000000' and a[2] == '0.957200':
					btype = int(a[0])

				# Check for TIP4P angle
				if a[1] == '7777.000000' and a[2] == '104.520000':
					atype = int(a[0])

			# Check if all types have been found. Return values if found
			if otype > -1 and htype > -1 and btype > -1 and atype > -1:
				return otype, htype, btype, atype


def PotEngSurfaceJob(run_name, input_script, system, domain, spanMolecule, resolution = .1,
					 queue='batch', procs=1, email='', pair_coeffs_included=True,
					 hybrid_pair=False,split=0,floor=[0.0,0.0,0.0]):
	"""
	Runs a set of LAMMPS jobs on the queue, writing one main data file and then
	modifying it incrementally to move one particular atom, in order to create
	a potential energy surface. This is faster than running write_lammps_data
	multiple times, since nearly all information is kept the same. Other than this,
	this function should be used almost exactly as job().
	This will dump all results into the lammps folder as separate folders because
	doing something else would require changing many other functions.
	
	Preconditions:
	run_name is a string of positive length, 31 or smaller.
	system is a valid utils.System instance.
	
	domain is a 3-element list of ints or floats, referring to the XYZ domain
	in which the molecule spanMolecule must span over. E.g., [0.0,1.0,1.0] will
	span over a domain which is 1x1 Angstrom box in the positive y and z directions
	from spanMolecule's original position.
	
	spanMolecule is a string which represents a path to a valid cml file, which
	represents the molecule which will be spanned across a surface in system.
	
	increment is an int or float which represents the resolution, in Angstroms,
	with which the atom will be spanned across the given domain. This span is
	linear.
	
	hybrid_pair is a boolean flag which should be set to True if the system
	should exhibit hybrid pairing interactions. (e.g. lj/cut and nm/cut)
	
	"""
	if len(run_name) > 31 and queue is not None:
		raise Exception("Job name too long (%d) for NBS. Max character length is 31." % len(run_name))

	# Change to correct directory
	os.system('mkdir -p lammps/%s' % run_name)
	os.chdir('lammps/%s' % run_name)
	
	#Generate the list of positions for the molecule to be tried
	vecList = _GenerateVecList(domain,resolution)
	
	#Find the position in the input script which declares the data file name
	readDataIndex = input_script.find("read_data")
	if readDataIndex == -1:
		raise ValueError("Invalid Input Script: no data read.")
	dotDataIndex = input_script.find(".data")
	
	
	#If the molecule is in the system, in its original position, remove it.
	if system.Contains(spanMolecule):
		system.Remove(spanMolecule)
	spanMolecule.SetCenter(floor)

	sys.setrecursionlimit(2000)
	#Dump to pickle files all global objects
	os.system('mkdir -p lammps/nbs_scripts')
	os.chdir('lammps/nbs_scripts')
	pickle.dump(system,open("system.pickle","wb"))
	pickle.dump(input_script,open("inputscript.pickle","wb"))
	pickle.dump(readDataIndex,open("readdataindex.pickle","wb"))
	pickle.dump(dotDataIndex,open("dotdataindex.pickle","wb"))
	pickle.dump(spanMolecule,open("spanmolecule.pickle","wb"))
	pickle.dump(pair_coeffs_included,open("pair_coeffs_included.pickle","wb"))
	pickle.dump(hybrid_pair,open("hybrid_pair.pickle","wb"))

	
	#For each and every vector in vecList, create a valid lammps job with the spanMolecule
	#traslated by that vector and run it. Between steps, remove the older spanMolecule
	#and add a new spanMolecule to the system which is in a new position, translated
	#by vec.
	vecList = [vecList]
	for vecLine in vecList:
		newVecs = []
		if split > 1:
			splitNum=len(vecLine)/split
			a = 0
			while a+splitNum < len(vecLine)-1:
				newVecs.append(vecLine[a:a+splitNum])
				a += splitNum
			newVecs.append(vecLine[a:])
		else:
			newVecs.append(vecLine)
		for b in range(len(newVecs)):
			lineName = run_name+"_"+str(vecLine[0][1])+"_sp"+str(b)
			_PotEngLine(lineName,run_name,input_script,spanMolecule,newVecs[b],procs=procs,
			email=email)
			#break
	


def PotEngSurfaceJob_v2(run_name, input_script, system, domain, spanMolecule, resolution = .1,
					 procs=1, email='', pair_coeffs_included=True,
					 hybrid_pair=False,assortBy='y',floor=[0.0,0.0,0.0]):
	"""
	Runs a set of LAMMPS jobs on the queue, writing one main data file and then
	modifying it incrementally to move one particular atom, in order to create
	a potential energy surface. This is faster than running write_lammps_data
	multiple times, since nearly all information is kept the same. Other than this,
	this function should be used almost exactly as job().
	
	This is a more optimized version of the original function, will run on a limited
	number of cores on the queue, run faster, and have lower overall CPU time.
	Rather than submitting all jobs independently to the queue, a single row of
	the PES in dimension assortBy is run locally on each processor. Thus, ~6400
	jobs running for 10 seconds each becomes ~80 jobs running for less than 800
	seconds each. 
	
	Preconditions:
	run_name is a string of positive length, 31 or smaller.
	system is a valid utils.System instance.
	
	domain is a 3-element list of ints or floats, referring to the XYZ domain
	in which the molecule spanMolecule must span over. E.g., [0.0,1.0,1.0] will
	span over a domain which is 1x1 Angstrom box in the positive y and z directions
	from spanMolecule's original position.
	
	spanMolecule is a string which represents a path to a valid cml file, which
	represents the molecule which will be spanned across a surface in system.
	
	increment is an int or float which represents the resolution, in Angstroms,
	with which the atom will be spanned across the given domain. This span is
	linear.
	
	hybrid_pair is a boolean flag which should be set to True if the system
	should exhibit hybrid pairing interactions. (e.g. lj/cut and nm/cut)
	
	assortBy is a string, one of ("x","y", or "z").
	
	"""
	if len(run_name) > 31 and queue is not None:
		raise Exception("Job name too long (%d) for NBS. Max character length is 31." % len(run_name))

	# Change to correct directory
	os.system('mkdir -p lammps/%s' % run_name)
	os.chdir('lammps/%s' % run_name)
	
	#Generate the list of positions for the molecule to be tried
	vecList = _GenerateVecList(domain, resolution, plotting=True, assortBy=assortBy)
	#Find the position in the input script which declares the data file name
	readDataIndex = input_script.find("read_data")
	if readDataIndex == -1:
		raise ValueError("Invalid Input Script: no data read.")
	dotDataIndex = input_script.find(".data")
	
	
	#If the molecule is in the system, in its original position, remove it.
	if system.Contains(spanMolecule):
		system.Remove(spanMolecule)
	spanMolecule.SetCenter(floor)
	
	sys.setrecursionlimit(2000)
	#Dump to pickle files all global objects
	os.system('mkdir -p lammps/nbs_scripts')
	os.chdir('lammps/nbs_scripts')
	pickle.dump(system,open("system.pickle","wb"))
	pickle.dump(input_script,open("inputscript.pickle","wb"))
	pickle.dump(readDataIndex,open("readdataindex.pickle","wb"))
	pickle.dump(dotDataIndex,open("dotdataindex.pickle","wb"))
	pickle.dump(spanMolecule,open("spanmolecule.pickle","wb"))
	pickle.dump(pair_coeffs_included,open("pair_coeffs_included.pickle","wb"))
	pickle.dump(hybrid_pair,open("hybrid_pair.pickle","wb"))
	
	#For every vector line in vecList, run a job to acquire a processor to run it,
	#and run a lammps job for each one.
	for vecLine in vecList:
		lineName = run_name+"_"+str(vecLine[0][1])
		_PotEngLine(lineName,run_name,input_script,spanMolecule,vecLine,procs=procs,
			email=email,pair_coeffs_included=pair_coeffs_included,hybrid_pair=hybrid_pair)
		#break
		


def _PotEngLine(lineName,run_name,input_script, spanMolecule, vecLine, procs=1, email='',
			pair_coeffs_included=True,hybrid_pair=False):
	"""
	A helper function to PotEngSurfaceJob_v2(). Writes a python script to be
	run by the ICSE servers and submits a valid NBS job to run that python script.
	The ICSE processor which is grabbed will then run a series of lammps runs,
	each of which involving translating spanMolecule by a vector in vecLine, and
	then resetting its position after the run.
	"""
	currdir = os.getcwd()
	pickle.dump(vecLine,open("vec_"+lineName+".pickle","wb"))
	pickle.dump(spanMolecule,open("spanMolecule_"+lineName+".pickle","wb"))
	
	NBS = open(lineName+".nbs",'w')
	NBStext='''#!/bin/sh
##NBS-name: \''''+lineName+'''\'
##NBS-nproc: 1
##NBS-queue: batch
##NBS-email: ''
source /fs/home/aec253/.zshrc

# Shared LAMMPS Library Configuration
export PYTHONPATH=/fs/home/yma3/Software/lammps_git/python:$PYTHONPATH
export PYTHONPATH=/fs/home/yma3/Projects/Forcefields/OPLS:$PYTHONPATH
export LD_LIBRARY_PATH=/fs/home/yma3/Software/lammps_git/src:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/mpich2/icse/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/fs/home/yma3/usr/local/lib/fftw3/lib:$LD_LIBRARY_PATH

/fs/home/aec253/anaconda/bin/python -u '''+lineName+'''.py >> '''+lineName+'''.log 2>&1
'''
	NBS.write(NBStext)
	NBS.close()
	
	pyScript = open(lineName+'.py','w')
	pyScriptText = '''from lammps_job import job
import cPickle as pickle
import os
from merlin import *
os.chdir(\''''+currdir+'''\')

run_name = "'''+run_name+'''"
system = pickle.load(open("system.pickle","rb"))
vecLine = pickle.load(open("vec_"+"'''+lineName+'''.pickle","rb"))
input_script = pickle.load(open("inputscript.pickle","rb"))
readDataIndex = pickle.load(open("readdataindex.pickle","rb"))
dotDataIndex = pickle.load(open("dotdataindex.pickle","rb"))
spanMolecule = pickle.load(open("spanmolecule.pickle","rb"))
pair_coeffs_included = pickle.load(open("pair_coeffs_included.pickle","rb"))
hybrid_pair = pickle.load(open("hybrid_pair.pickle","rb"))


os.chdir("..")
os.chdir("..")
spanLog = open("lammps/nbs_scripts/spanLog.txt","w")
for vec in vecLine:
	spanLog.write(str(vec)+"\\n")
	newRunName = run_name + "_" + str(vec[0])+ "x" + str(vec[1]) + "x" + str(vec[2])
	system.name = newRunName
	spanLog.write("Pre-Translated at "+`vec`+": "+`spanMolecule`+"\\n")
	newInputScript = input_script[:readDataIndex+10]+newRunName+input_script[dotDataIndex:]
	spanMolecule.translate(vec)
	spanLog.write("Post-Translated at "+`vec`+": "+`spanMolecule`+"\\n")
	system.add(spanMolecule)
	print "Running job " + newRunName
	job(newRunName, newInputScript, system, 
		pair_coeffs_included=pair_coeffs_included, hybrid_pair=hybrid_pair,
		hybrid_angle=False)
	system.Remove(spanMolecule)
	spanMolecule.translate([-vec[0],-vec[1],-vec[2]])
	#if vecLine.index(vec)>5:
		#break
spanLog.close()'''

	
	pyScript.write(pyScriptText)
	pyScript.close()
	
	os.system("jsub "+lineName+".nbs")


def OptJob(run_name, input_script, system, domain, spanMolecule, resolution = .1,
					 procs=1, email='', pair_coeffs_included=True,
					 hybrid_pair=False,assortBy='y',orientations=10,split=1,
					 floor=[0.0,0.0,0.0]):
	"""
	Runs a series of lammps jobs in a space defined by domain and resolution,
	at each position performing a rand_rotation of the spanMolecule. Runs the
	simulation at that position orienations many times. Each simulation line
	will be split over split many jobs. Other than this, this function can be
	used as PotEngSurfJob() and job().
	"""
	if len(run_name) > 31 and queue is not None:
		raise Exception("Job name too long (%d) for NBS. Max character length is 31." % len(run_name))
	
	# Change to correct directory
	os.system('mkdir -p lammps/%s' % run_name)
	os.chdir('lammps/%s' % run_name)
	
	#Generate the list of positions for the molecule to be tried
	vecList = _GenerateVecList(domain, resolution, plotting=False, assortBy=assortBy)
	#Find the position in the input script which declares the data file name
	readDataIndex = input_script.find("read_data")
	if readDataIndex == -1:
		raise ValueError("Invalid Input Script: no data read.")
	dotDataIndex = input_script.find(".data")
	
	
	#If the molecule is in the system, in its original position, remove it.
	if system.Contains(spanMolecule):
		system.Remove(spanMolecule)
	
	spanMolecule.SetCenter(floor)
	
	sys.setrecursionlimit(2000)
	#Dump to pickle files all global objects
	os.system('mkdir -p lammps/nbs_scripts')
	os.chdir('lammps/nbs_scripts')
	pickle.dump(system,open("system.pickle","wb"))
	pickle.dump(input_script,open("inputscript.pickle","wb"))
	pickle.dump(readDataIndex,open("readdataindex.pickle","wb"))
	pickle.dump(dotDataIndex,open("dotdataindex.pickle","wb"))
	pickle.dump(spanMolecule,open("spanmolecule.pickle","wb"))
	pickle.dump(pair_coeffs_included,open("pair_coeffs_included.pickle","wb"))
	pickle.dump(hybrid_pair,open("hybrid_pair.pickle","wb"))
	pickle.dump(orientations,open("orientations.pickle","wb"))
	
	#For every vector line in vecList, run a job to acquire a processor to run it,
	#and run a lammps job for each one.
	
	vecLine = vecList
	newVecs = []
	if split > 1:
		splitNum=len(vecLine)/split
		a = 0
		while a+splitNum < len(vecLine)-1:
			newVecs.append(vecLine[a:a+splitNum])
			a += splitNum
		newVecs.append(vecLine[a:])
	else:
		newVecs.append(vecLine)
	for b in range(len(newVecs)):
		lineName = run_name+"_"+str(vecLine[0][1])+"_sp"+str(b)
		_OptEngLine(lineName,run_name,input_script,spanMolecule,newVecs[b],procs=procs,
		email=email)
		#break


def _OptEngLine(lineName,run_name,input_script, spanMolecule, vecLine, procs=1, email='',
			orientations=10,split=1):
	"""
	A helper function to PotEngSurfaceJob_v2(). Writes a python script to be
	run by the ICSE servers and submits a valid NBS job to run that python script.
	The ICSE processor which is grabbed will then run a series of lammps runs,
	each of which involving translating spanMolecule by a vector in vecLine, and
	then resetting its position after the run.
	"""
	print "Running OptEngLine"
	currdir = os.getcwd()
	pickle.dump(vecLine,open("vec_"+lineName+".pickle","wb"))
	pickle.dump(spanMolecule,open("spanMolecule_"+lineName+".pickle","wb"))
	
	NBS = open(lineName+".nbs",'w')
	NBStext='''#!/bin/sh
##NBS-name: \''''+lineName+'''\'
##NBS-nproc: '''+str(procs)+'''
##NBS-queue: batch
##NBS-email: ''
source /fs/home/aec253/.zshrc

# Shared LAMMPS Library Configuration
export PYTHONPATH=/fs/home/yma3/Software/lammps_git/python:$PYTHONPATH
export PYTHONPATH=/fs/home/yma3/Projects/Forcefields/OPLS:$PYTHONPATH
export LD_LIBRARY_PATH=/fs/home/yma3/Software/lammps_git/src:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/mpich2/icse/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/fs/home/yma3/usr/local/lib/fftw3/lib:$LD_LIBRARY_PATH

/fs/home/aec253/anaconda/bin/python -u '''+lineName+'''.py >> '''+lineName+'''.log 2>&1
'''
	NBS.write(NBStext)
	NBS.close()
	
	pyScript = open(lineName+'.py','w')
	pyScriptText = '''from lammps_job import job
import cPickle as pickle
import os
from merlin import *
os.chdir(\''''+currdir+'''\')

run_name = "'''+run_name+'''"
system = pickle.load(open("system.pickle","rb"))
vecLine = pickle.load(open("vec_"+"'''+lineName+'''.pickle","rb"))
input_script = pickle.load(open("inputscript.pickle","rb"))
readDataIndex = pickle.load(open("readdataindex.pickle","rb"))
dotDataIndex = pickle.load(open("dotdataindex.pickle","rb"))
spanMolecule = pickle.load(open("spanmolecule.pickle","rb"))
pair_coeffs_included = pickle.load(open("pair_coeffs_included.pickle","rb"))
hybrid_pair = pickle.load(open("hybrid_pair.pickle","rb"))
orientations = pickle.load(open("orientations.pickle","rb"))


os.chdir("..")
os.chdir("..")
spanLog = open("lammps/nbs_scripts/spanLog.txt","w")
for vec in vecLine:
	for a in range(orientations):
		spanLog.write(str(vec)+"at orientation "+str(a)+"\\n")
		newRunName = run_name + "_" + str(vec[0])+ "x" + str(vec[1]) + "x" + (str(vec[2])
			+"_"+str(a))
		system.name = newRunName
		spanLog.write("Pre-Translated at "+`vec`+": "+`spanMolecule`+"\\n")
		spanMolecule.randRotateInPlace()
		newInputScript = input_script[:readDataIndex+10]+newRunName+input_script[dotDataIndex:]
		spanMolecule.translate(vec)
		spanLog.write("Post-Translated at "+`vec`+": "+`spanMolecule`+"\\n")
		system.add(spanMolecule)
		print "Running job " + newRunName
		job(newRunName, newInputScript, system, 
			pair_coeffs_included=pair_coeffs_included, hybrid_pair=hybrid_pair,
			hybrid_angle=False)
		system.Remove(spanMolecule)
		spanMolecule.translate([-vec[0],-vec[1],-vec[2]])
		#if vecLine.index(vec)>5:
			#break
spanLog.close()'''

	
	pyScript.write(pyScriptText)
	pyScript.close()
	
	os.system("jsub "+lineName+".nbs")


def _GenerateVecList(domain,resolution,plotting = False,assortBy='y'):
	"""
	A helper-function to generate a 2D-list of all translate vectors which should
	be tried by PotEngSurfaceJob. domain is an domain over which to be
	spanned (square prismatically) and resolution is the resolution at which the
	domain should be spanned.
	If plotting a 2d surface using this method, set plotting to True and the result
	will be a 3d-list, the outer of whihch will hold individual column's data.
	Precondition:
	domain must be a 3-element list of ints or floats,
	and resolution must be a positive int or float.
	plotting is a bool, or is not passed.
	"""
	#Initialize position lists
	xList= [0.0]
	yList= [0.0]
	zList= [0.0]
	vecList = []

	#Fill position lists with multiples of resolution up to the maximum domain size
	while abs(xList[-1]+resolution) <= abs(domain[0]):
		xList.append(xList[-1]+resolution)
	while abs(yList[-1]+resolution) <= abs(domain[1]):
		yList.append(yList[-1]+resolution)
	while abs(zList[-1]+resolution) <= abs(domain[2]):
		zList.append(zList[-1]+resolution)
	if not plotting:
		#Loop through all elements of all lists and add all of these points to posList
		for x in xList:
			for y in yList:
				for z in zList:
					vecList.append([x,y,z])
	elif assortBy=='y':
		for x in xList:
			for y in yList:
				addList = []
				for z in zList:
					addList.append([x,y,z])
				vecList.append(addList)
	elif assortBy == 'x':
		for z in zList:
			for x in xList:
				addList = []
				for y in yList:
					addList.append([x,y,z])
				vecList.append(addList)
	elif assortBy == 'z':
		for y in yList:
			for z in zList:
				addList = []
				for x in xList:
					addList.append([x,y,z])
				vecList.append(addList)
	
	return vecList





# A function to extract thermo output from lammps log files
# Automatically removes duplicate timesteps
# Adapted from log2txt.py from the lammps distribution
# Syntax:  log_file output_file X Y ...
#          log_file = LAMMPS log file
#          output_file = text file to create
#          X Y ... = columns to include (optional), X,Y are thermo keywords
#                    if no columns listed, all columns are included
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
# Modified by Yaset Acevedo to use updated log.py included in clancelot
def thermo_2_text(run_name, *properties):
	# Format log_file as needed
	log_file = 'lammps/%s/%s.log' % (run_name, run_name)
	if not os.path.isfile(log_file):
		raise IOError("No log file %s exists in %s." % (log_file, os.getcwd()))

	# Format output file name
	output_file = 'lammps/%s/thermo.txt' % (run_name)

	# Import log file
	lg = lammps_log(log_file)

	# If no properties specified, print out all properties
	if properties == []:
		lg.write(output_file)

	# Only print out selected properties
	else:
		str = "lg.write(output_file,"
		for word in properties:
			str += '"' + word + '",'
		str = str[:-1] + ')'
		eval(str)

# A function to return thermo output from lammps log files
def read_thermo(run_name, *properties):
	# Format log_file as needed
	log_file = 'lammps/%s/%s.log' % (run_name, run_name)
	if not os.path.isfile(log_file):
		raise IOError("No log file %s exists in %s." % (log_file, os.getcwd()))

	# Import log file
	lg = lammps_log(log_file)
	lg.last_modified = files.last_modified(log_path)

	# Print all property names
	names = lg.names
	txt = ''
	for name in names:
		txt += '%s ' % (name)
	print(txt)

	return lg
