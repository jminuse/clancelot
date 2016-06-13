import os, sys
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from utils import strip_colour, spaced_print

# Change font size
matplotlib.rcParams.update({'font.size': 10})

def is_int(x):
	try:
		int(x)
		return True
	except: return False

OLD_GRAPHS = None
#OLD_GRAPHS = './pys_chk/5/'
NEW_GRAPHS = './pys/'
LOOSE = 1
TO_GRAPH_OPTS = ['BFGS','LBFGS','SD','QM','FIRE']
TO_GRAPH_GEOM = ['CNH','CNLi','CNNa','CNK',
				 'BOH','BOLi','BONa','BOK']
SPACE_FIXING = {"LBFGS":"","BFGS":"  ","SD":"      ","QM":"     ","FIRE":"    "}

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
	output_file_1 = open("figs/%s_1" % geom,"w")
	output_file_2 = open("figs/%s_2" % geom,"w")
	out_F_2_data = []
	out_F_2_str = ""
	for opt in TO_GRAPH_OPTS:
		name = geom+"_"+opt
		if name in best_data_old:
			out_F_2_str += "%s_old\t" % opt
			min_y_list.append(min(best_data_old[name][1]))
			max_y_list.append(max(best_data_old[name][1]))
			x1, y1, l1 = range(len(best_data_old[name][1])), best_data_old[name][1], opt+"_old    step %d,    %.4f Ha/Ang,    %.2f kT_300" % (best_data_old[name][0], best_data_old[name][2], max_y_list[-1])
			out_F_2_data.append(y1)
			plt.plot(x1, y1, label=l1)
			output_file_1.write("%s\t%d\t%.4f\t%.2f\n" % (opt+"_old",best_data_old[name][0], best_data_old[name][2], max_y_list[-1]) )
		if name in best_data_new:
			min_y_list.append(min(best_data_new[name][1]))
			max_y_list.append(max(best_data_new[name][1]))
			if OLD_GRAPHS is None: legend_name_identifier = ""
			else: legend_name_identifier = "_new"
			out_F_2_str += "%s%s\t" % (opt,legend_name_identifier)
			x2, y2, l2 = range(len(best_data_new[name][1])), best_data_new[name][1], opt+"%s%s    step %d,    %.4f Ha/Ang,    %.2f kT_300" % (legend_name_identifier, SPACE_FIXING[opt], best_data_new[name][0], best_data_new[name][2], max_y_list[-1])
			plt.plot(x2, y2, label=l2)
			out_F_2_data.append(y2)
			output_file_1.write("%s\t%d\t%.4f\t%.2f\n" % (opt+legend_name_identifier,best_data_new[name][0], best_data_new[name][2], max_y_list[-1]) )
	output_file_1.close()

	out_F_2_str += "\n"
	for i in range(len(out_F_2_data[0])):
		for j in range(len(out_F_2_data)):
			out_F_2_str += "%.4f\t" % out_F_2_data[j][i]
		out_F_2_str += "\n"

	out_F_2_str = spaced_print(out_F_2_str,['\t'],buf=4)
	output_file_2.write(out_F_2_str)
	output_file_2.close()

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

	plt.savefig("./figs/%s" % (geom), dpi=100, bbox_inches='tight')
	plt.clf()
