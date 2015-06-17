import sys,  os
import g09, files
#import math, random, re,  utils, shutil, copy, cPickle

# From within a folder with a 'gaussian', scang will show energy landscape over
# several configurations, given logfiles 'name%.log' % step. using the command:
# scang ['name'|'gaussian/name'] [OPTIONAL: start step] [end step] [OPTIONAL: -t='title'] [OPTIONAL: -a="lowy highy"]

# ex. 'scang gaussian/ranthisjobnumer3- 2 7 -t="Energy Plot"'

name = sys.argv[1]
if name.startswith('gaussian/'):
	name=name[9:]

#custom axis?
if (sys.argv[len(sys.argv)-1]).startswith('-a='):
	axis=(sys.argv[len(sys.argv)-1])[3:]+'-'
	a=axis.split()
	print a
	ylow=int(a[0])
	yhigh=int((a[1])[:-1])
	sys.argv=sys.argv[:-1]
	custom_axis=True
else:
	# Default
	custom_axis=False

# look for title
if (sys.argv[len(sys.argv)-1]).startswith('-t='):
	title=(sys.argv[len(sys.argv)-1])[3:]
	sys.argv=sys.argv[:-1]
else:
	# Default plot title is the folder name
	title=os.path.basename(os.getcwd())

if len(sys.argv)==3:
	low=0
	count = int(sys.argv[2])
if len(sys.argv)==4:
	low = int(sys.argv[2])
	count = int(sys.argv[3])

f = open('out.xyz', 'w')
energies = []

print 'Step', 'E (Har)', 'Converged?'
for step in range(low,count):
	energy, atoms = g09.parse_atoms(name+str(step) if name.find('%')==-1 else name % step, check_convergence=False)
	files.write_xyz(atoms, f)
	print step, energy, int(g09.parse_atoms(name+str(step) if name.find('%')==-1 else name % step)!=None)
	energies.append(energy)

def matplot(y):
	import matplotlib.pyplot as plt
	plt.plot(y,marker='.')
	plt.xlabel('Step')
	plt.ylabel('E (kcal/mol)')
	plt.title(title)
	if custom_axis:
		plt.axis([low,count,ylow,yhigh])
	plt.show()

energies = [(e-energies[0])*627.5 for e in energies]
matplot(energies)
