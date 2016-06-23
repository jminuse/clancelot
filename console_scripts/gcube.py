from merlin import *
from subprocess import Popen

old_job = sys.argv[2]

if not os.path.exists('gaussian/%s.chk' % old_job):
	print 'Fatal error: file "gaussian/%s.chk" does not exist.' % old_job
	exit()

if not g09.parse_atoms(old_job):
	print 'Fatal error: "%s" is not converged. gcube does not work on unconverged jobs.' % old_job
	exit()

# Get the file to check
if not os.path.exists('gaussian/%s.fchk' % old_job):
	print 'Making gaussian/%s.fchk' % old_job
	Popen('/usr/local/gaussian/g09/g09/formchk gaussian/%s.chk gaussian/%s.fchk' % (old_job,old_job), shell=True).wait()

# Make the density and potential cube files
if not os.path.exists('gaussian/%s.cube' % (old_job+'_d')):
	print 'Making gaussian/%s.cube' % (old_job+'_d')
	Popen('/usr/local/gaussian/g09/g09/cubegen 0 density gaussian/%s.fchk gaussian/%s.cube 0 h'% (old_job,old_job+'_d'), shell=True).wait()

if not os.path.exists('gaussian/%s.cube' % (old_job+'_p')):
	print 'Making gaussian/%s.cube' % (old_job+'_p')
	Popen('/usr/local/gaussian/g09/g09/cubegen 0 potential gaussian/%s.fchk gaussian/%s.cube 0 h'% (old_job,old_job+'_p'), shell=True).wait()

if not (os.path.exists('gaussian/%s.cube' % (old_job+'_p')) and os.path.exists('gaussian/%s.cube' % (old_job+'_d')) ):
	print 'Fatal error: cube files not created'
	exit()

vmd_file = '''# Type logfile console into console to see all commands

# Get data
mol new gaussian/$$FPTR$$_d.cube
mol addfile gaussian/$$FPTR$$_p.cube

# Adjust first rep
mol modcolor 0 0 element
mol modstyle 0 0 CPK

# Adjust second rep
mol addrep 0
mol modcolor 1 0 Volume 1
mol modstyle 1 0 Isosurface 0.040000 0 0 0 1 1
mol modmaterial 1 0 Transparent'''.replace('$$FPTR$$',old_job)

f = open('tmp.vmd','w')
f.write(vmd_file)
f.close()

Popen('/fs/europa/g_pc/vmd-1.9 -e tmp.vmd', shell=True)
