import os, sys, subprocess
import files, g09, log

input = sys.argv[1]
try:
	vmd_open = sys.argv[2] in ['1','true','yes','y']
	print('sys.argv[2] = '+sys.argv[2])
	print('vmd_open = '+str(vmd_open))
except:
	vmd_open = True

if not input.endswith('.log'):
	input = 'gaussian/'+input+'.log'

output = sys.argv[3] if len(sys.argv)>3 else 'out'

energies, frames, time = g09.parse_atoms(input,parse_all=True)
for e in energies:
	print e

job_name = input[ input.rindex('/')+1 : input.rindex('.') ]

if time:
	print 'Job finished in %.2g seconds' % time
elif (job_name) in log.get_jlist():
	print 'Job is still running'
else:
	print 'Job failed to converge. Log file says:'
	os.system('tail -n 5 '+input)

if len(frames) > 0:
	files.write_xyz(frames, output)
	if vmd_open: os.system('/fs/europa/g_pc/vmd-1.9 '+output+'.xyz > /dev/null')
else:
	print 'There is no data in the log file'
