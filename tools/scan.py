import sys,  os, re
import g09, files
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



ex: 'scang gaussian/ranthisjobnumer3- 2 7 -t "Energy Plot" -x=0.000001,1.5423,3,...,7.0003452'
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
print 'Step', 'E (Har)', 'Converged?'
for step in range(low,count+1):
	energy, atoms = g09.parse_atoms(name+str(step) if name.find('%')==-1 else name % step, check_convergence=False)
	files.write_xyz(atoms, f)
	print step, energy, int(g09.parse_atoms(name+str(step) if name.find('%')==-1 else name % step)!=None)
	energies.append(energy)

def matplot(y):
	import matplotlib.pyplot as plt
	if not custom_x:
		plt.plot(y,marker='.')
	else:
		plt.plot(xs,y,marker='.')
	plt.xlabel(lx)
	plt.ylabel('E (kcal/mol)')
	plt.title(title)
	if custom_y:
		if not custom_x:
			plt.axis([low,count,ylow,yhigh])
		else:
			plt.axis([0,max(xs),ylow,yhigh])
	plt.show()

energies = [(e-energies[0])*627.5 for e in energies]
print energies
matplot(energies)
