import os, sys
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from utils import strip_colour

# Change font size
matplotlib.rcParams.update({'font.size': 10})

def is_int(x):
	try:
		int(x)
		return True
	except: return False

OLD_GRAPHS = './pys_chk/5/'
NEW_GRAPHS = './pys/'
LOOSE = 1
TO_GRAPH_OPTS = ['BFGS','LBFGS','SD','QM','FIRE']
TO_GRAPH_GEOM = ['CNH','CNLi','CNNa','CNK',
				 'BOH','BOLi','BONa','BOK']

best_data_new = {} # Will be stored as "CN$X_$OPT" or "BO$X_$OPT"
best_data_old = {}

if not os.path.isdir("./figs"): os.mkdir("figs")

# Read in best case data points of the new data set
for fptr in os.listdir(NEW_GRAPHS):
	LOOSE_FLAG = False
	if not fptr.endswith(".log"): continue
	file_by_line = open(NEW_GRAPHS+fptr,'r').read().split('\n')
	while len(file_by_line) > 0 and file_by_line[-1] == '': file_by_line = file_by_line[:-1]

	last_line = file_by_line[-1].split()
	while last_line == [] or not is_int(last_line[0]):
		file_by_line = file_by_line[:-1]
		if len(file_by_line) == 0:
			if LOOSE:
				print("Warning - No data for new graphs %s%s" % (NEW_GRAPHS, fptr))
				LOOSE_FLAG = True
				break
			else:
				raise Exception("Unable to find good last line in %s%s" % (NEW_GRAPHS, fptr))
		last_line = file_by_line[-1].split()
	if LOOSE_FLAG: continue

	# Parse last line for this data point
	data_name = fptr.split('_')[0] + '_' + fptr.split('.')[0].split('_')[-1]
	energy_data = [0.0]
	rms = None

	try:
		step = int(last_line[0])
		energy_data += [float(x) for x in last_line[3:-2]]
		rms = float(strip_colour(last_line[-2]))

		best_data_new[data_name] = [step, energy_data, rms]
	except:
		raise Exception("Unable to parse the last line in %s%s" % (NEW_GRAPHS, fptr))

# If it exists, read in best case data points from the previous run
if OLD_GRAPHS is not None:
	for fptr in os.listdir(OLD_GRAPHS):
		LOOSE_FLAG = False
		if not fptr.endswith(".log"): continue
		file_by_line = open(OLD_GRAPHS+fptr,'r').read().split('\n')
		while len(file_by_line) > 0 and file_by_line[-1] == '': file_by_line = file_by_line[:-1]

		last_line = file_by_line[-1].split()
		while last_line == [] or not is_int(last_line[0]):
			file_by_line = file_by_line[:-1]
			if len(file_by_line) == 0:
				if LOOSE:
					print("Warning - No data for old graphs %s%s" % (OLD_GRAPHS, fptr))
					LOOSE_FLAG = True
					break
				else:
					raise Exception("Unable to find good last line in %s%s" % (OLD_GRAPHS, fptr))

			last_line = file_by_line[-1].split()
		if LOOSE_FLAG: continue
		# Parse last line for this data point
		data_name = fptr.split('_')[0] + '_' + fptr.split('.')[0].split('_')[-1]
		energy_data = [0.0]
		rms = None

		try:
			step = int(last_line[0])
			energy_data += [float(x) for x in last_line[3:-2]]
			rms = float(strip_colour(last_line[-2]))

			best_data_old[data_name] = [step, energy_data, rms]
		except:
			raise Exception("Unable to parse the last line in %s%s" % (OLD_GRAPHS, fptr))

for geom in TO_GRAPH_GEOM:
	min_y_list = []
	max_y_list = []
	for opt in TO_GRAPH_OPTS:
		name = geom+"_"+opt
		if name in best_data_old:
			min_y_list.append(min(best_data_old[name][1]))
			max_y_list.append(max(best_data_old[name][1]))
			x1, y1, l1 = range(len(best_data_old[name][1])), best_data_old[name][1], opt+"_old    step %d,    %.4f Ha/Ang,    %.2f kT_300" % (best_data_old[name][0], best_data_old[name][2], max_y_list[-1])
			plt.plot(x1, y1, label=l1)
		if name in best_data_new:
			min_y_list.append(min(best_data_new[name][1]))
			max_y_list.append(max(best_data_new[name][1]))
			x2, y2, l2 = range(len(best_data_new[name][1])), best_data_new[name][1], opt+"_new    step %d,    %.4f Ha/Ang,    %.2f kT_300" % (best_data_new[name][0], best_data_new[name][2], max_y_list[-1])
			plt.plot(x2, y2, label=l2)

	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

	# Scale y-range to view everything
	axes = plt.gca()
	# Here we find the smallest one and add a buffer of 1.05 times it
	ymax = min(max_y_list) * 1.05
	ymin = min(min_y_list) * 1.05
	#ymax = sum(max_y_list)/float(len(max_y_list))
	axes.set_ylim([ymin,ymax])

	plt.xlabel("Reaction Coordinate")
	plt.ylabel("Energy Difference (kT_300)")
	plt.title("Isomerization of %s" % geom)
	fig = plt.gcf()

	#plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')

	plt.savefig("./figs/%s" % (geom), dpi=100, bbox_inches='tight')
	plt.clf()


'''

plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
plt.show()



t = np.arange(0.0, 2.0, 0.01)
s = np.sin(2*np.pi*t)
plt.plot(t, s)

plt.xlabel('time (s)')
plt.ylabel('voltage (mV)')
plt.title('About as simple as it gets, folks')
plt.grid(True)
plt.savefig("test.png")
plt.show()
'''