from merlin import *
from subprocess import Popen

old_job = sys.argv[1]

# Get the file to check
Popen('/usr/local/gaussian/g09/g09/formchk gaussian/%s.chk gaussian/%s.fchk' % (old_job,old_job), shell=True).wait()

# Make the density and potential cube files
Popen('/usr/local/gaussian/g09/g09/cubegen 0 density gaussian/%s.fchk gaussian/%s.cube 0 h'% (old_job,old_job+'_d'), shell=True).wait()
Popen('/usr/local/gaussian/g09/g09/cubegen 0 potential gaussian/%s.fchk gaussian/%s.cube 0 h'% (old_job,old_job+'_p'), shell=True).wait()

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
mol modmaterial 1 0 Transparent'''.replace('$$FPTR$$',sys.argv[1])

f = open('tmp.vmd','w')
f.write(vmd_file)
f.close()

Popen('/fs/europa/g_pc/vmd-1.9 -e tmp.vmd', shell=True)