import sys, os

run_name = sys.argv[1]
data = open(run_name+'.data').read()

start = data.index('Masses')
try:
	end = data.index('Pair Coeffs')
except:
	end = data.index('Bond Coeffs')

elements_by_mass = {12:'C', 14:'N', 1:'H', 16:'O', 32:'S', 35:'Cl', 207:'Pb'}
elements_by_index = {}

for line in data[start:end].splitlines():
	if line and line[0].isdigit():
		index, mass = line.split()
		elements_by_index[index] = elements_by_mass[int(round(float(mass)))]

f = open('out.xyz', 'w')
for line in open(run_name+'.xyz'):
	columns = line.split()
	if len(columns)>3:
		index, x, y, z = columns
		f.write("%s\t%s\t%s\t%s\n" % (elements_by_index[index], x, y, z))
	else: f.write(line)
f.close()

os.system('/fs/europa/g_pc/vmd-1.9 out.xyz > /dev/null')

