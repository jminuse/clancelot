import sys
sys.path.append('/fs/home/hch54/clancelot/tools/')
import files, neb

fptr = 'BOLi_LiBO.xyz'
frames = files.read_xyz('/fs/home/hch54/Documents/projects/BFGS_bench4/xyz/'+fptr)
opt = 'QM'
route = '! RI-B2PLYP D3BJ def2-TZVP def2-TZVP/C Grid3 FinalGrid5'

run_name = fptr[:fptr.find('.xyz')] + '_' + opt
neb.neb(run_name, frames, route, opt=opt, maxiter=1000, gtol=0.00183726383934, DFT='orca', alpha=0.1, dt=0.1, mem=40)