import sys
import files, neb

fptr = 'CNLi_LiCN'
frames = files.read_xyz('/fs/home/jms875/clancelot/unit_tests/neb/xyz/'+fptr+'.xyz')
opt = 'LBFGS'
route = '! HF-3c Grid3 FinalGrid5'

run_name = fptr[:fptr.find('.xyz')] + '_' + opt
neb.neb(run_name, frames, route, opt=opt, maxiter=200, gtol=0.00183726383934, DFT='orca', alpha=0.1, dt=0.1, mem=40, Nmax=20)