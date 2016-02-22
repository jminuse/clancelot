from merlin import *
from getpass import getuser
#############################################################################################
## FOR EASY USE - Choose a test bench using the commandline
## Goes from fast (1) to slow (5).
## NOTE!!! IF YOU WANT YOUR OWN PERSONAL TESTBENCH, ADD AT 6 ONWARDS!!!
use_testbench = 1

user = getuser()
os.chdir('/fs/home/%s/clancelot/unit_tests/neb/' % user)
if len(sys.argv) > 1: use_testbench = sys.argv[1]
queue = 'batch'

print("Using testbench %s" % use_testbench)

if use_testbench == '1':
	maxiter, gtol = 10, units.convert('eV/Ang','Ha/Ang',0.1)
	route = '! M062X def2-TZVP Grid3 FinalGrid5'
	opts = ['SD','BFGS','LBFGS']
	xyzs = ['CNH_HCN']

	## Ensure you maintain similarities between optimization methods
	alpha = 0.1
	dt = 0.1
	DFT = 'orca'
	mem = 40
	Nmax = 20
elif use_testbench == '2':
	maxiter, gtol = 30, units.convert('eV/Ang','Ha/Ang',0.1)
	route = '! M062X def2-TZVP Grid3 FinalGrid5'
	opts = ['QM','SD','FIRE','BFGS','LBFGS']
	xyzs = ['CNH_HCN','CNLi_LiCN','CNNa_NaCN','CNK_KCN']

	## Ensure you maintain similarities between optimization methods
	alpha = 0.1
	dt = 0.1
	DFT = 'orca'
	mem = 40
	Nmax = 20
elif use_testbench == '3':
	maxiter, gtol = 100, units.convert('eV/Ang','Ha/Ang',0.1)
	route = '! M062X def2-TZVP Grid3 FinalGrid5'
	opts = ['QM','SD','FIRE','BFGS','LBFGS']
	xyzs = ['CNH_HCN','CNLi_LiCN','CNNa_NaCN','CNK_KCN','BOH_HBO','BOLi_LiBO','BONa_NaBO','BOK_KBO']

	## Ensure you maintain similarities between optimization methods
	alpha = 0.1
	dt = 0.1
	DFT = 'orca'
	mem = 40
	Nmax = 20
elif use_testbench == '4':
	maxiter, gtol = 1000, units.convert('eV/Ang','Ha/Ang',0.05)
	route = '! M062X def2-TZVP Grid3 FinalGrid5'
	opts = ['QM','SD','FIRE','BFGS','LBFGS']
	xyzs = ['CNH_HCN','CNLi_LiCN','CNNa_NaCN','CNK_KCN','BOH_HBO','BOLi_LiBO','BONa_NaBO','BOK_KBO']

	## Ensure you maintain similarities between optimization methods
	queue = 'long'
	alpha = 0.1
	dt = 0.1
	DFT = 'orca'
	mem = 40
	Nmax = 20
elif use_testbench == '5':
	maxiter, gtol = 1000, units.convert('eV/Ang','Ha/Ang',0.05)
	route = '! RI-B2PLYP D3BJ def2-TZVP def2-TZVP/C Grid3 FinalGrid5'
	opts = ['QM','SD','FIRE','BFGS','LBFGS']
	xyzs = ['CNH_HCN','CNLi_LiCN','CNNa_NaCN','CNK_KCN','BOH_HBO','BOLi_LiBO','BONa_NaBO','BOK_KBO']

	## Ensure you maintain similarities between optimization methods
	queue = 'long'
	alpha = 0.1
	dt = 0.1
	DFT = 'orca'
	mem = 40
	Nmax = 20
else:
	print("Warning, choose a test bench available in clancelot/unit_tests/neb/run.py")
	sys.exit()

#############################################################################################

# A String for the simulation job
s_hold = '''#!/bin/bash
##NBS-name: "$$$$$$"
##NBS-nproc: 2
##NBS-queue: "$QUEUE$"

source /fs/home/$USER/.zshrc
/fs/home/$USER/anaconda/bin/python2.7 -u /fs/home/$USER/clancelot/unit_tests/neb/pys/$$$$$$.py > /fs/home/$USER/clancelot/unit_tests/neb/pys/$$$$$$.log 2>&1
'''

# For each NEB, throw on a queue
for fptr in xyzs:
	for opt in opts:
		s = '''import sys
import files, neb

fptr = '$FPTR$'
frames = files.read_xyz('/fs/home/$USER$/clancelot/unit_tests/neb/xyz/'+fptr+'.xyz')
opt = '$OPT$'
route = '$ROUTE$'

run_name = fptr[:fptr.find('.xyz')] + '_' + opt
neb.neb(run_name, frames, route, opt=opt, maxiter=$MAXITER$, gtol=$GTOL$, DFT='$DFT$', alpha=$ALPHA$, dt=$DT$, mem=$MEM$, Nmax=$NMAX$)'''

		# Replace with defined test
		s = s.replace('$FPTR$',fptr)
		s = s.replace('$OPT$',opt)
		s = s.replace('$ROUTE$',route)
		s = s.replace('$GTOL$',str(gtol))
		s = s.replace('$MAXITER$',str(maxiter))
		s = s.replace('$ALPHA$',str(alpha))
		s = s.replace('$DT$',str(dt))
		s = s.replace('$MEM$',str(mem))
		s = s.replace('$NMAX$',str(Nmax))
		s = s.replace('$USER$',user)
		s = s.replace('$DFT$',DFT)

		# Write the python file
		f = open('pys/'+fptr + '_' + opt+'.py','w')
		f.write(s)
		f.close()

		# Write the NBS file
		s = fptr + '_' + opt+'.py'
		f=open('pysub.nbs','w')
		f.write(s_hold.replace('$$$$$$',s[:-3]).replace("$QUEUE$",queue))
		f.close()

		# Run the simulation
		os.system('jsub pysub.nbs')

os.system('rm pysub.nbs')
