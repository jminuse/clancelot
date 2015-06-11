import os, sys
from fnmatch import fnmatch
import log

# Get jlist
jlist = ' '.join(log.get_jlist())

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