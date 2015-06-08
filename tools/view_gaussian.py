import os, sys, subprocess
import files, g09

input = sys.argv[1]
if not input.endswith('.log'):
	input = 'gaussian/'+input+'.log'

output = sys.argv[2] if len(sys.argv)>2 else 'out'

energies, frames, time = g09.parse_all(input)
for e in energies:
	print e

jlist = subprocess.Popen('jlist', shell=True, stdout=subprocess.PIPE).communicate()[0]
job_name = input[ input.rindex('/')+1 : input.rindex('.') ]

if time:
	print 'Job finished in %.2g seconds' % time
elif (' '+job_name+' ') in jlist:
	print 'Job is still running'
else:
	print 'Job failed to converge. Log file says:'
	os.system('tail -n 5 '+input)

if len(frames) > 0:
	files.write_xyz(frames, output)
	os.system('/fs/europa/g_pc/vmd-1.9 '+output+'.xyz > /dev/null')
else:
	print 'There is no data in the log file'
