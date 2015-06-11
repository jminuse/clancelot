import os, sys, re
import units, g09
from fnmatch import fnmatch
from subprocess import Popen, PIPE
from getpass import getuser

def get_jlist():
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

	return jlist.split()

# A function that checks if the gaussian job 'run_name' exists in gaussian.log
def chk_gaussian(run_name,sptr=None,force=False):
	if sptr == None: 
		try: sptr = open('gaussian.log').readlines()
		except: return -1

	# Loop through string
	for i,s in enumerate(sptr):
		s = s.split()
		# Find a header for a run
		if (len(s) == 1) and (s[0][-1]==':'):
			# Check if it's the same
			if(run_name == s[0][:-1]):
				# Print a warning if force == True, else throw error
				if force:
					#print('\nWarning! You already have a run called %s in your log file. Will append to run, but the gaussian files will be overwritten.\n' % run_name)
					return i
				else:
					raise Exception('\nYou already have a run called %s in your log file. If you want to re-write the gaussian data, please use force=True.\n' % run_name)
	return -1

# A function that puts information in the log file
def put_gaussian(run_name,route,extra_section,blurb,eRec,force=False):
	if blurb: blurb = blurb.replace('\n\n','\n')
	if route: route = route.replace('\n\n','\n---Empty Line---\n')
	chk = True
	try: s_hold = open('gaussian.log').readlines()
	except:
		f = open('gaussian.log','w')
		f.write('Simulation Name:\n\tGaussian Job Type\n\tDescription of the Simulation\n\tFinal System Energy\n')
		f.close()
		s_hold = open('gaussian.log').readlines()
		chk = False

	log_info = '\t'+route+'\n'
	try: log_info += '\t'+blurb+'\n'
	except: pass
	if extra_section != '':
		extra_section = extra_section.replace('\n\n','\n---Empty Line---\n')
		extra_section = extra_section.replace('\n','\n\t')
		log_info += '\t'+extra_section+'\n'
		log_info = log_info.replace('\n\t\n','\n')
	if eRec: log_info += '\t##$$@@#'+run_name+'#$@$#@#$\n'
	if log_info: log_info = log_info.replace('\n\n','\n')
	log_info += '\n'

	# If we want to append to a pre-existing file
	if chk:
		# We ensure we aren't re-writing a data file
		index = chk_gaussian(run_name,sptr=s_hold,force=force)
		# If we are re-writing and we are forcing it to do so, we append to the log
		if (index > -1):
			while index < len(s_hold):
				if s_hold[index] == '\n':
					break
				index += 1

		if (index > -1): s_hold[index] = '-------\n'+log_info
		if (index == -1) and (s_hold[-1] == '\n'): s_hold.append(run_name+':\n'+log_info)
		elif (index == -1): s_hold.append('\n'+run_name+':\n'+log_info)

	f = open('gaussian.log','w')
	f.write(''.join(s_hold))
	f.close()

# A function to get energy data for output to screen and/or record to gaussian.log
def chkg_E(fptr,unit='Ha',record=False,e_list=False,suppress=False):
	# Read in data from file
	energies, _, time = g09.parse_all("gaussian/"+fptr+".log")

	# If you want the standard output to terminal, do this
	if not suppress:
		# Get all energy values
		if e_list:
			for e in energies: print(e)
		print('\n---------------------------------------------------')
		print('Energy Data Points: '+str(len(energies)))
		if len(energies)>2: print('dE 2nd last = '+str(units.convert_energy('Ha',unit,energies[-2]-energies[-3]))+' '+unit)
		if len(energies)>1: print('dE last = '+str(units.convert_energy('Ha',unit,energies[-1]-energies[-2]))+' '+unit)
		if len(energies)>0: print('Last Energy = '+str(energies[-1])+' Ha')
		print('---------------------------------------------------')
		if time: print 'Job finished in %.2g seconds' % time
		elif (' '+fptr+' ') in get_jlist(): print 'Job is still running'
		else: 
			print 'Job failed to converge. Log file says:\n~~~~ End Of File Info'
			os.system('tail -n 5 '+"gaussian/"+fptr+".log")
			print '~~~~ Convergenge Criteria'
			s = open('gaussian/'+fptr+'.log').read()
			print('\n'.join(s[s.rfind("Converged?"):].split('\n')[1:5]))
		print('---------------------------------------------------\n')
	# If you want to record data, do this
	if record:
		try: s = open('gaussian.log').read().replace('##$$@@#'+fptr+'#$@$#@#$',str(units.convert_energy('Ha',unit,energies[-1]))+' '+unit)
		except:
			s = open('gaussian.log').read().replace('##$$@@#'+fptr+'#$@$#@#$','Error - Could not get the energy from file')
			print('\nWarning - Could not get the energy for '+fptr+'.\n')
		
		f = open('gaussian.log','w')
		f.write(s)
		f.close()

# Goes through and updates all the eval values that haven't been updated
def chkg_E_all(u='Ha'):
	sptr = open('gaussian.log').read()
	a = '##$$@@#'
	b = '#$@$#@#$'
	s = sptr[sptr.find(a)+len(a):sptr.find(b)]
	if (sptr.find(a)==-1) or (sptr.find(b)==-1): return 0
	while len(s)>0:
		chkg_E(s,record=True,unit=u,suppress=True)
		sptr = open('gaussian.log').read()
		if (sptr.find(a)==-1) or (sptr.find(b)==-1): return 0
		s = sptr[sptr.find(a)+len(a):sptr.find(b)]
