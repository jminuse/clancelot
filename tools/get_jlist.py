import os, sys
from fnmatch import fnmatch
from subprocess import Popen, PIPE
from getpass import getuser

USER_NAME=getuser()

# Get input from jlist as a string
p = Popen(['jlist'], stdout=PIPE)
output = p.stdout.read().split()

# Make an empty string to hold our list
jlist = ""

# Loop through jlist and get file names
for i,s in enumerate(output):
	if s==USER_NAME:
		jlist = jlist + output[i+1] + " "

# Take care of using wildcards in searching jlist.
# For instance, if I have 'test1' 'test2' and 'ethane' running and want to list 'test*'
if len(sys.argv) > 1:
		qlist = ''
		for s in jlist.split():
			if(fnmatch(s,sys.argv[-1])):
				qlist = qlist + s + " "
		jlist = qlist

# Return your list for script processing
print jlist