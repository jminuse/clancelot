import os, sys
import files, g09

input = sys.argv[1]
if not input.endswith('.log'):
	input = 'gaussian/'+input+'.log'

output = sys.argv[2] if len(sys.argv)>2 else 'out'

energies, frames, time = g09.parse_all(input)
for e in energies:
	print e

if time:
	print 'Time = %.2g seconds' % time
else:
	print 'Job is not converged. Log file says:'
	os.system('tail -n 5 '+input)

files.write_xyz(frames, output)
os.system('/fs/europa/g_pc/vmd-1.9 '+output+'.xyz > /dev/null')

