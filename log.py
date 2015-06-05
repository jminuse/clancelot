import os, sys, re
import units

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
					print('\nWarning! You already have a run by this name in your log file. Will append to run, but the gaussian files will be overwritten.\n')
					return i
				else:
					raise Exception('\nError! You have this run already in the data file. If you want to re-write the gaussian data, please use force=True.\n')
	return -1

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

def chkg_E(fptr,record=False,unit='Ha',suppress=False):
	contents = open('gaussian/'+fptr+'.log').read()

	if not suppress:
		print '\n'.join(re.findall('SCF Done: +\S+ += +(\S+)', contents))
		print '\n'
		print('Final Energy = '+str(units.convert_energy('Ha',unit,float(re.findall('SCF Done: +\S+ += +(\S+)', contents)[-1])))+' '+unit+'\n')

		if 'Normal termination of Gaussian 09' not in contents:
			print 'Job did not finish'
		else:
			m = re.search('Job cpu time: +(\S+) +days +(\S+) +hours +(\S+) +minutes +(\S+) +seconds', contents)
			time = float(m.group(1))*24*60*60 + float(m.group(2))*60*60 + float(m.group(3))*60 + float(m.group(4))
			print m.group(0)
			print "%.2e s" % time

		if 'Counterpoise: corrected energy =' in contents:
			print 'Counterpoise E = ',
			print '\n'.join(re.findall('Counterpoise: corrected energy = +(\S+)', contents))

	if record:
		try: s = open('gaussian.log').read().replace('##$$@@#'+fptr+'#$@$#@#$',str(units.convert_energy('Ha',unit,float(re.findall('SCF Done: +\S+ += +(\S+)', contents)[-1])))+' '+unit)
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
