import os, sys, re
import units, g09
from fnmatch import fnmatch
from subprocess import Popen, PIPE
from getpass import getuser

def get_jlist(verbose=0):
	# Get input from jlist as a string
	p = Popen(['jlist'], stdout=PIPE)
	output = p.stdout.read()

	# Get data from string
	pattern = getuser()+'''[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)'''
	info = re.findall(pattern,output)

	# Get a list of names
	names = []
	for a in info: names.append(a[0])

	# If user wants more information
	if verbose == 2: 
		for i,a in enumerate(info):
			p = Popen(['jshow',a[0]], stdout=PIPE)
			s = p.stdout.read() 
			serv = s[s.find('Queue name:'):].split()[2].strip()
			try: threads = s[s.find('Slot Reservations'):].split()[4].strip()
			except: threads=1
			info[i] = info[i] + (serv,threads,)
		return info

	# Return appropriate information
	if verbose == 1: return info
	
	return names

# A function that checks if today's log file exists and if it does returns it as a string, if
# it doesn't it moves all old log files to the logs folder and generates a new file, returning it as
# a string. chk simply says False if we need to generate a new file.
def get_log_file():
	import time
	today_log = 'gaussian_%s.log' % time.strftime("%d_%m_%Y")

	chk = True
	try: s_hold = open(today_log).readlines()
	except:
		# Make sure log folder exists
		try:
			if not os.path.isdir('logs'): os.mkdir('logs')
		except: 
			print('\nWARNING - Could not make log folder...\n')

		# Move all older log files to log folder
		os.system('mv gaussian_*.log logs 2>/dev/null')

		# Open new log file for today
		f = open(today_log,'w')
		f.write('Simulation Name:\n\tGaussian Job Type\n\tDescription of the Simulation\n\tFinal System Energy\n')
		f.close()
		s_hold = open(today_log).readlines()
		chk = False

	return s_hold,chk

# A function that checks a string (sptr) for the logged data for run_name
def chk_log_str(sptr,run_name,force,fptr):
	# Throwing in a pre-check for older log files in gaussian directory in-case user is new to using logging.
	if os.path.exists('gaussian/'+run_name+'.log'):
		if force: pass
		else: 
			print('\nERROR - Log file gaussian/%s.log exists. If you want to re-write the gaussian data, please use force=True.\n' % run_name)
			sys.exit(-1)

	# Loop through string
	for i,s in enumerate(sptr):
		s = s.split()
		# Find a header for a run
		if (len(s) == 1) and (s[0][-1]==':'):
			# Check if it's the same
			if(run_name == s[0][:-1]):
				# Throw error if force != True, else do it
				if force:
					return i
				else:
					print('\nERROR - You already have a log entry for %s in %s. If you want to re-write the gaussian data, please use force=True.\n' % (run_name,fptr))
					sys.exit(-1)
	return i + 1

# A function that checks if the gaussian job 'run_name' exists in gaussian.log files
# Returns -2 for read errors
# Returns -1 if log data doesn't exist
# Returns 0 -> n for index of where the overlap is. Index is in regards to line
def chk_gaussian(run_name,sptr=None,force=False,neb=[False,None,None,None]):
	# String for today's log file
	import time
	today_log = 'gaussian_%s.log' % time.strftime("%d_%m_%Y")

	# Depending on neb get the right string to look for
	if neb[0]: log_name = neb[1]
	else: log_name = run_name

	# sptr holds data from today's log file (or is a string passed by user)
	if sptr == None: 
		try: sptr,_ = get_log_file()
		except: return -2,None

	# Check if string contains log_name
	index = chk_log_str(sptr,log_name,force,today_log)
	if index != len(sptr): return index,today_log

	# Now check all older log files
	for fptr in os.listdir('logs'):
		name, ext = os.path.splitext(fptr)
		if ext.lower() == '.log':
			fptr = os.path.join('logs/', fptr)
			try: sptr = open(fptr).readlines()
			except: return -2,None
			index = chk_log_str(sptr,log_name,force,fptr)
			if index != len(sptr): return index,fptr

	return -1,None

# A function that puts information in the log file
def put_gaussian(run_name,route,extra_section='',blurb=False,eRec=True,force=False,neb=[False,None,None,None]):
	# This section is parsing the data to look nice for output ----------------------------------------
	if blurb: blurb = blurb.replace('\n\n','\n')
	if route: route = route.replace('\n\n','\n---Empty Line---\n')

	# Note, neb data is: neb=[True/False, Name with %d, Number of frames, Current Frame]
	if neb[0]: eRec = False

	log_info = '\t'+route+'\n'
	try: log_info += '\t'+blurb+'\n'
	except: pass
	if extra_section != '':
		extra_section = extra_section.replace('\n\n','\n---Empty Line---\n')
		extra_section = extra_section.replace('\n','\n\t')
		log_info += '\t'+extra_section+'\n'
		log_info = log_info.replace('\n\t\n','\n')
	if eRec: log_info += '\t##$$@@#'+run_name+'#$@$#@#$\n'
	elif neb[0]: log_info += '\tNumber of frames = [0,'+str(neb[2]-1)+']\n'
	if log_info: log_info = log_info.replace('\n\n','\n')
	log_info += '\n'
	# -------------------------------------------------------------------------------------------------
	
	# String for today's log file
	import time
	today_log = 'gaussian_%s.log' % time.strftime("%d_%m_%Y")

	# Try opening the gaussian log file. If you can't, then make one
	s_hold,chk = get_log_file()

	w_chk = True
	# If neb, we only need to record once
	if neb[0] and (''.join(s_hold).find(neb[1])>-1) and neb[3]>0:
		chk = False
		w_chk = False

	# We want a string for the log file
	if neb[0]: log_name = neb[1]
	else: log_name = run_name

	# If we want to append to a pre-existing file
	if chk:
		# We ensure we aren't re-writing a data file
		index,fptr = chk_gaussian(log_name,sptr=s_hold,force=force,neb=neb)

		# ----------- Case 1: We have data in today's log file, append to that data
		if ((index > -1) and (fptr == today_log)):
			# Loop to the first newline
			while index < len(s_hold):
				if s_hold[index] == '\n': break
				index += 1
			s_hold[index] = '-------\n'+log_info
		# ----------- Case 2: We have data in older log file, just save data normally with fptr
		elif ((index > -1) and (fptr != today_log)):
			log_info += '\tolder_data in '+fptr+'\n'
			log_info = log_info.replace('\n\n','\n') + '\n'
			if s_hold[-1] == '\n': s_hold.append(log_name+':\n'+log_info)
			else: s_hold.append('\n'+log_name+':\n'+log_info)
		# ----------- Case 3: We have a new data set, just save data normally
		elif (index == -1):
			if s_hold[-1] == '\n': s_hold.append(log_name+':\n'+log_info)
			else: s_hold.append('\n'+log_name+':\n'+log_info)

	if w_chk:
		f = open(today_log,'w')
		f.write(''.join(s_hold))
		f.close()

# A function to get energy data for output to screen and/or record to gaussian.log
def chkg_E(fptr,unit='Ha',record=False,e_list=False,suppress=False):
	import time
	today_log = 'gaussian_%s.log' % time.strftime("%d_%m_%Y")
	# Read in data from file
	energies, _, time = g09.parse_atoms("gaussian/"+fptr+".log",parse_all=True)

	# If you want the standard output to terminal, do this
	if not suppress:
		# Get all energy values
		if e_list:
			for e in energies: print(e)
		print('\n---------------------------------------------------')
		print('Job Name: '+fptr)
		print('Energy Data Points: '+str(len(energies)))
		if len(energies)>2: print('dE 2nd last = '+str(units.convert_energy('Ha',unit,energies[-2]-energies[-3]))+' '+unit)
		if len(energies)>1: print('dE last = '+str(units.convert_energy('Ha',unit,energies[-1]-energies[-2]))+' '+unit)
		if len(energies)>0: print('Last Energy = '+str(energies[-1])+' Ha')
		print('---------------------------------------------------')
		if time: print 'Job finished in %.2g seconds' % time
		elif (fptr) in get_jlist(): 
			print 'Job is still running'
			print '~~~~ Convergenge Criteria'
			s = open('gaussian/'+fptr+'.log').read()
			print('\n'.join(s[s.rfind("Converged?"):].split('\n')[1:5]))
		else: 
			print 'Job failed to converge. Log file says:\n~~~~ End Of File Info'
			os.system('tail -n 5 '+"gaussian/"+fptr+".log")
			print '~~~~ Convergenge Criteria'
			s = open('gaussian/'+fptr+'.log').read()
			print('\n'.join(s[s.rfind("Converged?"):].split('\n')[1:5]))
		print('---------------------------------------------------\n')
	# If you want to record data, do this
	if record:
		try: s = open(today_log).read().replace('##$$@@#'+fptr+'#$@$#@#$',str(units.convert_energy('Ha',unit,energies[-1]))+' '+unit)
		except:
			s = open(today_log).read().replace('##$$@@#'+fptr+'#$@$#@#$','Error - Could not get the energy from file')
			print('\nWarning - Could not get the energy for '+fptr+'.\n')
		
		f = open(today_log,'w')
		f.write(s)
		f.close()

# Goes through and updates all the eval values that haven't been updated
def chkg_E_all(u='Ha'):
	import time
	today_log = 'gaussian_%s.log' % time.strftime("%d_%m_%Y")
	sptr = open(today_log).read()
	a = '##$$@@#'
	b = '#$@$#@#$'
	s = sptr[sptr.find(a)+len(a):sptr.find(b)]
	if (sptr.find(a)==-1) or (sptr.find(b)==-1): return 0
	while len(s)>0:
		chkg_E(s,record=True,unit=u,suppress=True)
		sptr = open(today_log).read()
		if (sptr.find(a)==-1) or (sptr.find(b)==-1): return 0
		s = sptr[sptr.find(a)+len(a):sptr.find(b)]
