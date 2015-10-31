import sys,  os, re
import g09, orca, files
from merlin import units

#import math, random,  utils, shutil, copy, cPickle
help = '''USE: scang name start end [OPTION]*

From within a folder with a 'gaussian', scang will show energy landscape over
several configurations, given logfiles 'name%.log' % step

INPUTS:
name     :   can either be 'name' or 'gaussian/name'
start    :   integer starting step (optional, DEFAULT=0)
end      :   integer end step

OPTIONS:
-h , --help                 : print this message and exit
-t 'title' , -t='title'     : sets the plot title to the custom 
                               title (DEFAULT=current directory name)
-lx 'x label', -lx='label'  : sets x-axis label (DEFAULT='Step')
-x fl,...,fl , -x=fl,...,fl : sets the x-coordinates for the n steps 
                               to the given values instead of graphing
                               energy vs. step
-y low,high , -y=low,high   : sets the two endpoints for the y axis to
                               the given values
-c char,start,stop          : given a character, replaces it from start
                               to stop and runs scang on each.  Graphs
                               are compiled together.
-dft type                   : given a string describing what to scan. By
                               default this is 'g09' but can be 'orca'.
-u units                    : given a unit for y axis. By default this
                               is kcal/mol
-neb value                  : force neb endpoint on (1) or off (2) default
                              is 0 for don't force anything.

ex: 'scang gaussian/ranthisjobnumer3- 2 7 -t "Energy Plot" -x=0.000001,1.5423,3,...,7.0003452'
ex: 'scang gaussian/test_neb_$-%d_yay 1 10 -c $,0,6'

Note, in the second example it will include 0 and 11 for the full range of
neb data.
'''
#Print help message if too few arguments
if len(sys.argv) <3:
	# print 'ERROR: not enough arguments\n\n'
	print help
	raise SystemExit

#accepts either the log name or gaussian/name
name = sys.argv[1]
if name.startswith('gaussian/'):
	name=name[9:]
if name.startswith('orca/'):
	name=name[5:]

#Clean up flags to make parsing easier
arg = []
for s in sys.argv:
	if s.startswith('--'):
		s=s[1:]
	if '=' in s:
		arg = arg + s.split('=')
	else:
		arg.append(s)

#Default Settings:
custom_x=False
custom_y=False
lx='Step'
title=os.path.basename(os.getcwd())
comp=None
compE=[]
s_units='kcal/mol'
neb_force=0
DFT='g09'

#Parse optional flags
unidentifiable=[]
if (len(arg)>=3):
	curr_length=None
	old_length=len(arg)
	while curr_length != old_length and len(arg) > 2:
		old_length=len(arg)

		#Help message:
		if arg[-1].startswith('-h'):
			print help
			raise SystemExit
		#Custom Title (default is the current directory name)
		if arg[-2].startswith('-t'):
			print "  Custom Title = " + arg[-1]
			title=arg[-1]
			arg=arg[:-2]
		#Custom y-axis endpoints (default is autoscale)
		elif arg[-2].startswith('-y'):
			print "  Custom y-axis endpoints = " + arg[-1]
			ys=arg[-1].split(',')
			for i in range(len(ys)):
				ys[i]=float(ys[i])
			custom_y=True
			if len(ys)!=2:
				print "\n\nERROR: custom y axis flag given, but two endpoints not given"
				raise SystemExit
			ylow=ys[0]
			yhigh=ys[1]
			arg=arg[:-2]
		#Custom x-coordinates (default is step number)
		elif arg[-2].startswith('-x'):
			print "  Custom x coordinates = " + arg[-1]
			xs=arg[-1]
			xs=xs.split(",")
			for i in range(len(xs)):
				xs[i]=float(xs[i])
			custom_x=True
			arg=arg[:-2]
		#Custom x-axis labels (default is 'Step')
		elif arg[-2].startswith('-lx'):
			print "  Custom x axis label = " + arg[-1]
			lx=arg[-1]
			arg=arg[:-2]
		#Custom range for compiling together data
		elif arg[-2].startswith('-c'):
			print " Custom compilation = " + arg[-1]
			comp=arg[-1].split(',')
			arg=arg[:-2]
		#Custom y axis units
		elif arg[-2].startswith('-u'):
			print " Units = " + arg[-1]
			s_units=arg[-1]
			arg=arg[:-2]
		#NEB compile
		elif arg[-2].startswith('-neb'):
			print " Neb Force = " + str(arg[-1])
			neb_force= int(arg[-1])
			arg=arg[:-2]
		#Scan over g09 or orca
		elif arg[-2].startswith('-dft'):
			print " DFT scan = " + arg[-1]
			DFT=arg[-1]
			if DFT not in ['g09','orca']:
				print("Error - DFT must be 'g09' or 'orca', not %s" % DFT)
				sys.exit()
			arg=arg[:-2]

		curr_length=len(arg)

		#If the flag wasn't identified, add to list
		if curr_length==old_length and curr_length > 4:
			unidentifiable.insert(0,arg[-1])
			curr_length=curr_length-1
			arg=arg[:-1]

	#Print unidentified flags
	if unidentifiable!=[]:
		print "\nCannot identify the following flag(s)/input(s), please refer to the help documentation via -h or --help"
		print "The following will be ignored: "+str(unidentifiable) + '\n'

# If we are compiling together the data, we want to loop between the different data sets
if comp != None: loops=range(int(comp[1]),int(comp[2])+1)
else: loops = [1]

for loop in loops:
	#Set low and count values, based on whether low was supplied or not
	if len(sys.argv)==3 or (len(sys.argv) > 3 and not re.match('\d+',sys.argv[3])):
		low=0
		count = int(sys.argv[2])
	if len(sys.argv)>=4 and re.match('\d+',sys.argv[3]):
		low = int(sys.argv[2])
		count = int(sys.argv[3])

	#If custom x-coords given, ensure that x and y values are supplied 1-1
	if custom_x and (len(range(low,count+1))) != len(xs):
		if (len(range(low,count))) > len(xs):
			print '\n\nERROR: applied \'-x\' flag with too few x-coordinates specified for the number of frames to be plotted.'
			print '%d frames to be plotted vs. %d x coordinates given' % ((len(range(low,count))) , len(xs))
			raise SystemExit
		else:
			print '\n\nERROR: applied \'-x\' flag with too many x-coordinates specified for the number of frames to be plotted.'
			print '%d frames to be plotted vs. %d x coordinates given' % ((len(range(low,count))) , len(xs))
			raise SystemExit
	f = open('out.xyz', 'w')
	energies = []

	#Print and parse energy values
	if comp == None: print 'Step', 'E (Har)', 'Converged?'

	# Here we want to take into account having to collect an extra data point for neb calculations
	t_low, t_count = low, count
	if ((comp != None and loop == 0) or (neb_force==1 and loop == loops[0])) and neb_force!=2:
		t_low -= 1
		t_count += 1

	# Here we will get the energy data
	if ((comp != None and loop > 0) or (neb_force==1 and loop == loops[0])) and neb_force!=2: energies.append(offset_low)
	for step in range(t_low,t_count+1):
		tmp_name = name+str(step) if name.find('%')==-1 else name % step # Get step
		if comp != None: tmp_name = tmp_name.replace(comp[0],str(loop)) # Get loop
		if DFT == 'g09':
			energy, atoms = g09.parse_atoms(tmp_name, check_convergence=False)
		elif DFT == 'orca':
			atoms, energy = orca.parse_atoms(tmp_name)
		files.write_xyz(atoms, f)
		if comp == None: 
			if DFT == 'g09':
				print step, energy, int(g09.parse_atoms(tmp_name)!=None)
			elif DFT == 'orca':
				print step, energy
		energies.append(energy)
		if ((comp != None and loop == 0) or (neb_force==1 and loop == loops[0])) and neb_force!=2:
			if step == t_low: offset_low = energy
			if step == t_count: offset_high = energy
	if ((comp != None and loop > 0) or (neb_force==1 and loop == loops[0])) and neb_force!=2: energies.append(offset_high)

	if comp == None:
		def matplot(y):
			import matplotlib.pyplot as plt
			if not custom_x:
				plt.plot(y,marker='.')
			else:
				plt.plot(xs,y,marker='.')
			plt.xlabel(lx)
			plt.ylabel('E (%s)' % s_units)
			plt.title(title)
			if custom_y:
				if not custom_x:
					plt.axis([t_low,t_count,ylow,yhigh])
				else:
					plt.axis([0,max(xs),ylow,yhigh])
			plt.show()

		energies = [units.convert_energy('Ha',s_units, e-energies[0]) for e in energies]
		print energies
		matplot(energies)
	else:
		if ((comp != None and loop > 0) or (neb_force==1 and loop == loops[0])) and neb_force!=2:
			energies = [units.convert_energy('Ha',s_units, e-offset_low) for e in energies]
		else:
			energies = [units.convert_energy('Ha',s_units, e-energies[0]) for e in energies]
		compE.append(energies)


def comp_matplot(yy,start_val):
	import matplotlib.pyplot as plt
	if not custom_x:
		for i,y in enumerate(yy):
			plt.plot(y,marker='.',label=str(int(start_val) + i))
	else:
		for y in yy:
			plt.plot(xs,y,marker='.',label=str(int(start_val) + i))
	plt.xlabel(lx)
	plt.ylabel('E (%s)' % s_units)
	plt.title(title)
	if custom_y:
		if not custom_x:
			plt.axis([low,count,ylow,yhigh])
		else:
			plt.axis([0,max(xs),ylow,yhigh])
	plt.legend()
	plt.show()

# Get max and min, and generate plot
if comp != None:
	a = float('-inf')
	b = float('inf')
	for c in compE:
		if max(c) > a: a = max(c)
		if min(c) < b: b = min(c)
	print "Max E = %lg, Min E = %lg (Units = %s)\n" % (a,b,s_units)

	comp_matplot(compE,comp[1])
