import os, sys
from fnmatch import fnmatch

s_ext = sys.argv[1]

pylist = ''
d = sys.argv[2]+'/'

for fptr in os.listdir(d):
	name, ext = os.path.splitext(fptr)
	if ext.lower() == s_ext:
		pylist += name+s_ext+' '

if len(sys.argv) > 3:
	for i,s in enumerate(sys.argv):
		if i<3: continue
		print(s+' ')
else: print pylist



#import os, sys
#from fnmatch import fnmatch
#pylist = ''
#d = sys.argv[1]+'/'
#for fptr in os.listdir(d):
#	name, ext = os.path.splitext(fptr)
#	if ext.lower() == '.py':
#		pylist += name+'.py '
#if len(sys.argv) > 2:
#	for i,s in enumerate(sys.argv):
#		if i<2: continue
#		print(s+' ')
#else: print pylist