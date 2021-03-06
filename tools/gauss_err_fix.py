# To fully install run the following command to open the crontab editor (vi):
#			crontab -e
# Now in here type the following (adjusting accordingly):
#			MAILTO=""
#			10 * * * * /usr/bin/python /fs/home/hch54/clancelot/tools/gauss_err_fix.py

USER_PATH = ['/fs/home/hch54/Documents/projects/NanoCrystal/gaussian/']

import sys, os
sys.path.append('/fs/home/hch54/clancelot/tools')
os.environ["PATH"] += ':/opt/voyager/nbs/bin/:/usr/common/bin/'
from merlin import *
from ast import literal_eval as l_eval
import shutil
end_path = os.path.dirname(os.path.realpath(__file__))

def process_err(fptr,err,job_info):
	# This error arises when one gets 180 degree angles.  To fix this you need to resubmit the job with a slightly new geometry
	if err == 'FormBX had a problem.':
		# Get the routing line from a previous job
		route,_ = g09.parse_route(open(fptr+'.inp').readline()[2:].strip())

		# Get the atom list from previous job

		# If the previous run had '/gen', we need to copy that for an extra_section
		if(route[0].split('/')[1].lower() in ['gen','genecp']):
			extras = open(fptr+'.inp').read().split('\n\n')[3:]
			extras = '\n\n'.join(extras)
		else: extras = ''

		# Append geom new def if we need to
		need_geom = True
		for i,p in enumerate(route):
			if 'geom' in p.lower(): need_geom = False		
		if need_geom: route.append('GEOM=(Check,NewDefinition)')

		# Get previous job info
		serv = job_info[3]
		threads = job_info[4]

		# Make sure we won't re-write a job with new job name
		path = fptr[:fptr.rfind('/')+1]
		i = 0
		while os.path.exists(path+job_info[0]+'_R'+str(i)+'.log'): i+= 1
		job_name = job_info[0]+'_R'+str(i)
		
		# Copy chkpoint file
		shutil.copyfile(fptr+'.chk', path+job_name+'.chk')

		# Submit new job
		old_path = os.getcwd()
		path = path[:path[:-1].rfind('/')+1] # Move back another directory
		os.chdir(path)
		g09.job(job_name, ' '.join(route), queue=serv, extra_section=extras ,blurb='Re_running %s: %s' % (job_info[0],err),procs=threads, previous=job_info[0], force=True, err=True)

def get_err(fptr):
	flag = False
	for u in USER_PATH:
		try:
			path = u+fptr
			err = open(path+'.log').read()
			err = ''.join(err[:err.rfind('Error termination via Lnk1e')].split('\n')[-2:]).strip()
			flag = True
		except: pass
	if flag: return path,err
	else: return None,None # Return nothing if we couldn't find anything

# Get previous list of all runs
try: old=l_eval(open(end_path+'/jlist.list').read())
except: old=[]
# Get new list
new=log.get_jlist(2)
# Find finished jobs
for a in old:
	found = False
	for b in new:
		if a[0] == b[0]:
			found = True
	if not found:
		fptr,err = get_err(a[0])
		process_err(fptr,err,a)

# Re-write list for jobs running
f = open(end_path+'/jlist.list','w')
f.write(str(new))
f.close()