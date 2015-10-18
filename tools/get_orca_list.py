import os, sys
from fnmatch import fnmatch

gauss_list = ''
d = sys.argv[1]+'/orca/' # Get the directory you want to search

try:
	# Loop through and get a string list of all files of the desired extension
	for fptr in os.listdir(d):
		name = os.path.splitext(fptr)
		orca_list += name+' '

	# If you want to use wildcards (*), this takes care of that
	if len(sys.argv) > 2:
		flist = ''
		for s in orca_list.split():
			if(fnmatch(s,sys.argv[-1])):
				flist = flist + s + " "
		orca_list = flist

	# Return your list for script processing
	print orca_list
except: pass
