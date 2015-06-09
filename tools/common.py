import os, sys
import math

# These imports are for the notify and watch functions -----
from subprocess import Popen, PIPE
from time import sleep
import smtplib
from email.mime.text import MIMEText
# ----------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# These two functions help notify me when a simulation job has finished
def notify(job_name, cell):
	msg = MIMEText(job_name + " has finished")

	msg['Subject'] = job_name + " has finished"
	msg['From'] = "g_pc@cbe.cornell.edu"
	msg['To'] = cell

	# Send the message via our own SMTP server, but don't include the envelope header.
	s = smtplib.SMTP('localhost')
	s.sendmail("g_pc@cbe.cornell.edu", [cell], msg.as_string())
	s.quit()

def watch(job_name, cell):
	while 1:
		p = Popen(['jlist'], stdout=PIPE)
		s = p.stdout.read().split()
		dne = True
		for i in range(len(s)):
			if((s[i] == job_name) and ((i + 2) < len(s)) and (s[i+2] == 'RUNNING')):
				dne = False
		if(dne):
			break
		sleep(10)

	notify(job_name, cell)

# --------------------------------------------------------------------------------------------------------------
# This function will get the data for a time correlation so we can calculate the standard error
# t_break should be given in femtoseconds.  t_len is the length of a timestep and t_step is how many time steps per output
def stdErr(data, stats, t_len, t_step, t_break=-1, d_block=1, start=0, stop=None):
	# P = t_b * sigma_b ^ 2 / (<A^2> - <A>^2)
	P_denom = stats['stdev']**2
	p_inv = []
	t_inv = []
	if not stop: stop = len(data)

	# Loop through some block sizes
	for i in range(start, stop, d_block):
		# Loop through blocks
		n = int(math.floor(len(data)/i))
		variance_b = 0

		# Here we will sum up the variance for n blocks of size i
		for j in range(n):
			# Sum the data element-wise.  Output is now a list
			sums = [0.,0.]
			for k in range(i):
				sums = map(sum,zip(sums,data[j*i+k][0]))

			#variance_b = sums[1] - sums[0]**2
			avg_b = sums[0]/i
			variance_b += (avg_b - stats['avg'])**2

		# Finalize the variance_b calculation by dividing by the number of blocks
		variance_b /= n

		# Here we have t_len femtoseconds per each timestep and i*t_step timesteps per block thermo output
		t_inv.append(1./float(t_len*i*t_step))
		# Note, we are doing inverse so P = t_b * variance_b / P_denom becomes P_denom / (t_b * variance_b)
		p_inv.append(P_denom / ((t_len*i*t_step)*variance_b))

		# Stop doing things around 1 ns because roughly here onwards we'll get statistically insignificant data
		if ((t_break != -1) and (t_len*i*t_step > t_break)):
			break

	return t_inv, p_inv
# --------------------------------------------------------------------------------------------------------------
# Generate a histogram data set.  Pass data, bin size, start (optional), stop (optional)
def get_hist(data, step=None, num=None, start=None, stop=None, norm=False, buf=1):
	if not ((step==None) != (num==None)):
		print step
		print num
		raise Exception('Must input either bin step or number of bins.  Not neither, and not both')

	bins = [0]
	counts = [0]
	tot = 0

	# Get max and min vals
	max_v = float("-inf")
	min_v = float("inf")

	for i in data:
		if float(max_v) < i: max_v = i
		if float(min_v) > i: min_v = i
		if norm: tot += 1

	if not start: start = min_v
	if not stop: stop = max_v
	if not num: num = int((max_v - min_v)/step) + buf
	if not step: step = float((stop - start)/float(num))

	bins = [i*step+start for i in range(num)]
	count = [0]*len(bins)

	for d in data:
		count[int((d-min_v) / step)] += 1 # This is [ ) type histogram

	if norm:
		for i in range(len(count)):
			count[i] /= float(tot*step)

	return bins, count

# --------------------------------------------------------------------------------------------------------------
# Generate a 2D histogram data set.  Pass data, bin size
def get_hist2D(data, step=None, num=None, start=None, stop=None, norm=False, buf=1, cont=None, fptr=None, ftype=None):
	if not ((step==None) != (num==None)):
		print step
		print num
		raise Exception('Must input either bin step or number of bins.  Not neither, and not both')

	bins_x = []
	bins_y = []
	if not cont: count = [[0]]
	else: count = cont
	tot = 0.

	# Get max and min vals
	max_x = float("-inf")
	min_x = float("inf")
	max_y = float("-inf")
	min_y = float("inf")

	for i in data:
		if float(max_x) < i[0]: max_x = i[0]
		if float(min_x) > i[0]: min_x = i[0]
		if float(max_y) < i[1]: max_y = i[1]
		if float(min_y) > i[1]: min_y = i[1]
		if norm: 
			if not cont: tot += 1.

	if type(buf) is float:
		a = buf
		buf = [0]*2
		buf[0],buf[1] = a,a
	if type(buf) is int:
		a = buf
		buf = [0]*2
		buf[0],buf[1] = a,a

	if type(num) is float:
		a = num
		num = [0]*2
		num[0],num[1] = a,a
	if type(num) is int:
		a = num
		num = [0]*2
		num[0],num[1] = a,a

	if type(step) is float:
		a = step
		step = [0]*2
		step[0],step[1] = a,a
	if type(step) is int:
		a = step
		step = [0]*2
		step[0],step[1] = a,a

	if type(start) is float:
		a = start
		start = [0]*2
		start[0],start[1] = a,a
	if type(start) is int:
		a = start
		start = [0]*2
		start[0],start[1] = a,a

	if type(stop) is float:
		a = stop
		stop = [0]*2
		stop[0],stop[1] = a,a
	if type(stop) is int:
		a = stop
		stop = [0]*2
		stop[0],stop[1] = a,a

	if start != None:
		min_x = start[0]
		min_y = start[1]

	if stop != None:
		max_x = stop[0]
		max_y = stop[1]

	if not num: 
		num = [0]*2
		num[0] = int((max_x - min_x)/step[0]) + buf[0]
		num[1] = int((max_y - min_y)/step[0]) + buf[1]
	if not step: 
		step = [0]*2
		step[0] = float((max_x - min_x)/float(num[0]))
		step[1] = float((max_y - min_y)/float(num[1]))

	bins_x = [i*step[0]+min_x for i in range(num[0])]
	bins_y = [i*step[1]+min_y for i in range(num[1])]
	if not cont: count = [[0 for j in range(len(bins_y))] for i in range(len(bins_x))]

	for d in data:
		try:
			count[int((d[0]-min_x) / step[0])][int((d[1]-min_y) / step[1])] += 1 # This is [ ) type histogram
		except IndexError:
			print('Index Error has occured:\n')
			print('\td[0] = '+str(d[0])+'\n')
			print('\tmin_x = '+str(min_x)+'\n')
			print('\tstep[0] = '+str(step[0])+'\n')
			print('\td[1] = '+str(d[1])+'\n')
			print('\tmin_y = '+str(min_y)+'\n')
			print('\tstep[1] = '+str(step[1])+'\n')
			print('\tlen(count) = '+str(len(count))+'\n')
			print('\tlen(coun[0]) = '+str(len(count[0]))+'\n')
			print('\tlen(bins_x) = '+str(len(bins_x))+'\n')
			print('\tlen(bins_y) = '+str(len(bins_y))+'\n')
			print('\ttried x:  = '+str(int((d[0]-min_x) / step[0]))+'\n')
			print('\ttried y:  = '+str(int((d[1]-min_y) / step[1]))+'\n')
			return
	if norm:
		if cont!=None:
			tot = 0.
			for i in range(len(count)):
				for j in range(len(count[i])):
					tot += count[i][j]

		for i in range(len(count)):
			for j in range(len(count[i])):
				try:
					count[i][j] /= float(tot*step[0]*step[1])
				except ZeroDivisionError:
					print('Tried dividing by zero:\n')
					print('\ttot = '+str(tot)+'\n')
					print('\tstep[0] = '+str(step[0])+'\n')
					print('\tstep[1] = '+str(step[1])+'\n')
					return

	if fptr:
		f = open(fptr,'w')
		if not ftype:
			for row in count:
				for val in row:
					f.write(str(val)+'\t')
				f.write('\n')
		elif ftype == 'gnuplot':
			for i,x in enumerate(bins_x):
				for j,y in enumerate(bins_y):
					f.write(str(x)+'\t'+str(y)+'\t'+str(count[i][j])+'\n\n')
		elif ftype == 'matrix':
			f.write('NaN\t')
			for x in bins_x:
				f.write(str(x)+'\t')
			f.write('\n')
			for j,y in enumerate(bins_y):
				f.write(str(y)+'\t')
				for i,x in enumerate(bins_x):
					f.write(str(count[i][j])+'\t')
				f.write('\n')
		f.close()
	return zip(bins_x,bins_y), count

# --------------------------------------------------------------------------------------------------------------
# A function that gets all data from a column.  It takes a filename, an id string for a line before data, and what column to get data for
# Data is returned as a list of tuples holding (data, data^2), the total average, and the total stdev
def get_a_col_data(filename, s_id, c, t=1, bias=[[(None,None)]], end=False, start=0, stop=None, raw=False):

	#--------------------------------------------------------------------------
	# Ensure everything is in the right construction (mainly ensure bias is list of list of tuples)
	if type(c) is int:
		c = [c]
	# If only int/float, assume bias by multiplication and make into tuple
	done = False
	while done == False:
		done = True
		if type(bias) is float:
			done = False
			bias = (bias,'*')
		if type(bias) is int:
			done = False
			bias = (bias,'*')
		# If only tuple, make into list
		if type(bias) is tuple:
			done = False
			bias = [bias]
		# If only list, make into list of lists
		try:
			if type(bias[0]) is list:
				pass
			else:
				done = False
				bias = [bias]
		except:
			done = False
			bias = [bias]

		# If length is different than columns, fix that
		if (len(bias)!=len(c)):
			done = False
			for i in range(len(bias),len(c)):
				bias.append(bias[0])
			if(len(c) < len(bias)):
				bias = bias[:(len(c)-len(bias))]
	#--------------------------------------------------------------------------

	# Get the file
	f = open(filename,'r')

	# Loop through until the data
	for i in range(t):
		while True:
			try:
				while (f.readline().split()[0] != s_id): pass
			except EOFError: raise Exception('ERROR - Could not find the desired line')
			except: continue
			break

	# Initialize stuff
	data = []
	stats = [{'avg':None,'stdev':None,'max':None,'min':None,'sum':None,'len':None} for i in range(len(c))] # avg, stdev, max, min, sum, len
	avg = [0.]*len(c)
	stdev = [0.]*len(c)
	max_val = [(float("-inf"),0)]*len(c)
	min_val = [(float("inf"),0)]*len(c)
	line = -1

	# Loop through the data
	for l in f:

		line += 1
		if (line < start):
			continue

		try:
			if ((stop != None) and (line > stop)):
				continue
		except:
			pass

		col = l.split()

		# Initialize num_c once
		try:
			num_c
		except NameError:
			num_c = len(col)

		# Ensure we have data to parse
			# We have num_c lines
			# First value is numeric
		# Loop through the columns
		for j in range(len(c)):

			if((len(col)==num_c) and isNum(col[0])):
				if(j==0):
					
					# Deal with the bias
					x = float(col[c[j]])
					for b in bias[j]:
						if b[1]=='+':
							x += float(b[0])
						elif b[1]=='-':
							x -= float(b[0])
						elif b[1]=='/':
							x /= float(b[0])
						elif b[1]=='*':
							x *= float(b[0])
						elif b[1]=='%':
							x %= float(b[0])
					tmp_data = ((x, x**2),)

					if(tmp_data[0][0] > max_val[0][0]):
						max_val[0] = (tmp_data[0][0],line)
					if(tmp_data[0][0] < min_val[0][0]):
						min_val[0] = (tmp_data[0][0],line)

					avg[0] += tmp_data[0][0]
					stdev[0] += tmp_data[0][1]
				else:

					# Deal with the bias
					x = float(col[c[j]])
					for b in bias[j]:
						if b[1]=='+':
							x += float(b[0])
						elif b[1]=='-':
							x -= float(b[0])
						elif b[1]=='/':
							x /= float(b[0])
						elif b[1]=='*':
							x *= float(b[0])
						elif b[1]=='%':
							x %= float(b[0])
					tmp_data += ((x, x**2),)

					if(tmp_data[j][0] > max_val[j][0]):
						max_val[j] = (tmp_data[j][0],line)
					if(tmp_data[j][0] < min_val[j][0]):
						min_val[j] = (tmp_data[j][0],line)

					avg[j] += tmp_data[j][0]
					stdev[j] += tmp_data[j][1]

		if tmp_data == None: pass
		else: data.append(tmp_data)
		tmp_data = None

		if(end and len(col) != num_c):
			break

	# We use the fact that variance = stdev^2 = <x^2> - <x>^2 = avg(x^2) - avg(x)^2
	# That way we only need to loop through the data once
	# Note, if we expect data to sum up to a value larger than 10^308 (which I hope not), we will be approaching the limit of python's float
	# So we would need to do a rolling average/stdev in batches.  But this is highly unlikely right now so has not been implemented

	if(len(data) == 0):
		return None, None, False

	for i,stat in enumerate(stats):
		stat['sum'] = avg[i]
		stat['max'] = max_val[i]
		stat['min'] = min_val[i]
		stat['len'] = len(data)
		stat['avg'] = avg[i] / float(len(data))
		stat['stdev'] = math.sqrt(stdev[i]/float(len(data)) - (avg[i]/float(len(data)))**2)


	f.close()

	if raw:
		for i,d in enumerate(data):
			tmp = [0]*len(d)
			for j in range(len(d)):
				tmp[j] = d[j][0]
			data[i] = tuple(tmp)

	return data, stats, line > stop

# --------------------------------------------------------------------------------------------------------------
# Simple checks on numeric status of string
def isInt(s):
    try: 
        int(s)
        return True
    except ValueError: return False

def isNum(s):
    try: 
        float(s)
        return True
    except ValueError: return False

