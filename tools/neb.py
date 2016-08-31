from merlin import *
import re, shutil, copy

from scipy.linalg import block_diag
import numpy as np

# Clancelot Nudged Elastic Band package
# Currently supports g09
# Cite NEB: http://scitation.aip.org/content/aip/journal/jcp/113/22/10.1063/1.1323224
# BFGS is the best method, cite http://theory.cm.utexas.edu/henkelman/pubs/sheppard08_134106.pdf
# Nudged Elastic Band. k for VASP is 5 eV/Angstrom, ie 0.1837 Hartree/Angstrom.
# Gtol of 1E-5 from scipy.bfgs package
# Headers from http://patorjk.com/software/taag/#p=display&f=Banner3&t=RUN%20SIMULATION, font "Banner3"

##########################################################################
##########################################################################
##    ## ######## ########      ######   #######  ########  ######## 
###   ## ##       ##     ##    ##    ## ##     ## ##     ## ##       
####  ## ##       ##     ##    ##       ##     ## ##     ## ##       
## ## ## ######   ########     ##       ##     ## ##     ## ######   
##  #### ##       ##     ##    ##       ##     ## ##     ## ##       
##   ### ##       ##     ##    ##    ## ##     ## ##     ## ##       
##    ## ######## ########      ######   #######  ########  ########
##########################################################################
##########################################################################

FAIL_CONVERGENCE = -1
MAXITER_CONVERGENCE = 0
FMAX_CONVERGENCE = 1
GTOL_CONVERGENCE = 2


def g09_start_job(NEB, i, state, procs, queue, force, initial_guess, extra_section, mem):
    if NEB.step>0:
        guess = ' Guess=Read'
    else:
        if initial_guess:
            guess = ' Guess=Read'
        else:
            guess = '' #no previous guess for first step
    return g09.job('%s-%d-%d'%(NEB.name,NEB.step,i), NEB.theory+' Force'+guess, state, procs=procs, queue=queue, force=force, previous=('%s-%d-%d'%(NEB.name,NEB.step-1,i)) if NEB.step>0 else initial_guess, extra_section=extra_section+'\n\n', neb=[True,'%s-%%d-%%d'%(NEB.name),len(NEB.states),i], mem=mem)


def g09_results(NEB, step_to_use, i, state):
    result = g09.parse_atoms('%s-%d-%d' % (NEB.name, step_to_use, i), check_convergence=False, parse_all=False)
    if not result:
        raise Exception('parse_atoms failed')
    new_energy, new_atoms = result
    
    # Check if coordinates are aligned properly between state and new_atoms
    def check_atom_coords(atoms1,atoms2,precision=1e-6):
        for a1,a2 in zip(atoms1,atoms2):
            if abs(a1.x-a2.x)>precision or abs(a1.y-a2.y)>precision or abs(a1.z-a2.z)>precision:
                print i, 'atoms not in same frame:', a1.x, a1.y, a1.z, 'vs', a2.x, a2.y, a2.z
                print abs(a1.x-a2.x), abs(a1.y-a2.y), abs(a1.z-a2.z)
                exit()

    if i!=0 and i!=len(NEB.states)-1:
        check_atom_coords(state,new_atoms)
        for a,b in zip(state, new_atoms):
            a.fx = units.convert('Ha/Bohr','Ha/Ang',b.fx)
            a.fy = units.convert('Ha/Bohr','Ha/Ang',b.fy)
            a.fz = units.convert('Ha/Bohr','Ha/Ang',b.fz)    
    
    return new_energy, new_atoms


def orca_start_job(NEB, i, state, procs, queue, force, initial_guess, extra_section, mem):
    if NEB.step>0:
        previous = '%s-%d-%d' % (NEB.name,NEB.step-1,i)
    else:
        if initial_guess:
            if hasattr(initial_guess, '__iter__'):
                previous = initial_guess[i]
            else:
                previous = initial_guess
        else:
            previous = None
    return orca.job('%s-%d-%d'%(NEB.name,NEB.step,i), NEB.theory, state, extra_section=extra_section, grad=True, procs=procs, queue=queue, previous=previous)


def orca_results(NEB, step_to_use, i, state):
    new_atoms, new_energy = orca.engrad_read('%s-%d-%d' % (NEB.name, step_to_use, i),force='Ha/Ang',pos='Ang')
    
    for a,b in zip(state, new_atoms):
        a.fx,a.fy,a.fz = b.fx,b.fy,b.fz
    
    return new_energy, new_atoms


def neb(name, states, theory, extra_section='', spring_atoms=None, procs=1, queue=None,
        disp=0, k=0.1837, fit_rigid=True, start_job=None, get_results=None,
        DFT='orca', opt='BFGS', gtol=1e-3, fmax=None, maxiter=1000,
        alpha=0.1, beta=0.5, tau=1E-3, reset=10, H_reset=True, Nmax=20,
        viscosity=0.1, dtmax=1.0, Nmin=5, finc=1.1, fdec=0.5, astart=0.1, fa=0.99,
        step_min=1E-8, step_max=0.2, bt_max=None, linesearch='backtrack', L2norm=True, bt_eps=1E-3,
        dt = 0.3, euler=True, force=True, mem=25, blurb=None, initial_guess=None, start_from=0): 
    
    DFT = DFT.lower().strip()
    if DFT=='orca':
        start_job = orca_start_job
        get_results = orca_results
    if DFT=='g09':
        start_job = g09_start_job
        get_results = g09_results

    # Contemplating a force addition of NoSymm to route if (1) DFT='g09' and (2) NoSymm not said

    # Try seeing if neb was run for <= 2 frames
    if type(states) == list and type(states[0]) != list:
        print("Error - Only one frame in NEB calculation. Did you mean to run an optimization instead?")
        sys.exit()
    elif type(states) == type(states[0]) and len(states) <= 2:
        print("Error - NEB requires at least 3 frames to run. You have entered only %d frames." % len(states))
        sys.exit()   

    # Set which atoms will be affected by virtual springs
    if not spring_atoms: # If not given, select all
        spring_atoms = range(len(states[0]))
    elif type(spring_atoms)==str: # A list of element names
        elements = spring_atoms.split()
        spring_atoms = [i for i,a in enumerate(states[0]) if a.element in elements]

    # Output for user
    if opt == 'BROYDEN_ROOT':
        print("\nRunning neb with optimization method %s" % str(opt))   
    elif opt in ['QM','FIRE']:
        if euler: tmp = ' with euler'
        else: tmp = ' with verlet alg.'
        print("\nRunning neb with optimization method %s%s" % (str(opt), tmp))
        print("\tdt = %lg" % dt)
        if opt == 'QM':
            print("\tviscosity = %lg, step_max = %lg" % (viscosity,step_max))
        if opt == 'FIRE':
            print("\tdtmax = %lg, step_max = %lg" % (dtmax, step_max))
            print("\tNmin = %lg, finc = %lg, fdec = %lg" % (Nmin, finc, fdec))
            print("\tastart = %lg, fa = %lg" % (astart, fa))
    elif opt in ['BFGS','LBFGS']:
        print("\nRunning neb with optimization method %s" % str(opt))
        print("\talpha = %lg, beta = %lg" % (alpha, beta)),
        print("")
        print("\tH_reset = %s" % str(H_reset)),
        print(", reset = %s, linesearch = %s" % (str(reset), linesearch))
        print("\tstep_max = %lg" % step_max),
        print("")
    elif opt == 'SD':
        print("\nRunning neb with optimization method %s" % str(opt))
        print("\talpha = %lg" % alpha)
    else:
        print("\nERROR - %s optimizations method does not exist! Choose from the following:" % str(opt))
        print("\t1. BFGS")
        print("\t2. LBFGS")
        print("\t3. QM")
        print("\t4. SD")
        print("\t5. FIRE")
        print("\t6. BROYDEN_ROOT\n")
        sys.exit()

    print("Spring Constant for NEB: %lg Ha/Ang = %lg eV/Ang" % (k, units.convert_energy("Ha","eV",k)))
    print("Convergence Criteria: gtol = %lg (Ha/Ang)" % gtol)
    print("Run_Name = %s" % str(name))
    print("DFT Package = %s" % DFT)
    if blurb != None: print("-----\n" + blurb)
    print("---------------------------------------------")

    # Class to contain working variables
    class NEB:
        name, states, theory, k = None, None, None, None
        error, gradient = None, None
        step = 0
        def __init__(self, name, states, theory, extra_section='', k=0.1837, fmax=None, start_from=0):
            NEB.name = name
            NEB.states = states
            NEB.theory = theory
            NEB.extra_section = extra_section
            NEB.k = k
            NEB.prv_RMS = None
            NEB.prv_MAX = None
            NEB.prv_MAX_E = None
            NEB.nframes = len(states)
            NEB.RMS_force = float('inf')
            NEB.MAX_force = float('inf')
            NEB.MAX_energy = float('inf')
            NEB.FMAX = fmax
            NEB.step = start_from

            # In all cases, optimize the path via rigid rotations first
            # Set to only if fit_rigid for comparison purposes
            if fit_rigid: utils.procrustes(NEB.states) 
            
            # Load initial coordinates into flat array for optimizer
            NEB.coords_start = []
            for s in states[1:-1]:
                for a in s:
                    NEB.coords_start += [a.x, a.y, a.z]
    
        @staticmethod
        def calculate(coords):
            # Update coordinates in states.
            # This won't change anything on the first run through, but will on subsequent ones
            coord_count = 0
            for s in NEB.states[1:-1]:
                for a in s:
                    a.x, a.y, a.z = coords[coord_count], coords[coord_count+1], coords[coord_count+2]
                    coord_count += 3
            
            # Start DFT jobs
            running_jobs = []            
            for i,state in enumerate(NEB.states):
                if (i==0 or i==len(NEB.states)-1) and NEB.step>0:
                    pass # No need to calculate anything for first and last states after the first step
                else:
                    running_jobs.append( start_job(NEB, i, state, procs, queue, force, initial_guess, extra_section, mem) )
            # Wait for jobs to finish
            for j in running_jobs: j.wait()

            # Get forces and energies from DFT calculations
            energies = []
            for i,state in enumerate(NEB.states):
                # State 0 and state N-1 don't change, so just use result from NEB.step == 0
                if (i==0 or i==len(NEB.states)-1):
                    step_to_use = 0
                else:
                    step_to_use = NEB.step

                new_energy, new_atoms = get_results(NEB, step_to_use, i, state)
                energies.append(new_energy)

            # V = potential energy from DFT. energies = V+springs
            V = copy.deepcopy(energies)

            fake_states = copy.deepcopy(NEB.states)
            if fit_rigid: full_rotation = utils.procrustes(fake_states, append_in_loop=False)

            # DEBUGGING CODE HERE
            #f_chk = []
            #for state in NEB.states[1:-1]:
                #for a in state:
                    #f_chk += [a.fx, a.fy, a.fz]
            #for i,f in enumerate(f_chk): f_chk[i] = units.convert_energy("Ha","eV",f)
            #f_chk = np.array(f_chk).reshape((-1,3))
            #e_chk = [units.convert_energy("Ha","eV",e) for e in energies]
            #print "Check 1: Force, Energy = ", f_chk[3], ", ", e_chk[:5]
            
            # Get positions in a flat array
            def get_positions(image):
                pos = np.array([np.empty([3]) for j in image])
                for j,atom in enumerate(image):
                    if j not in spring_atoms: continue
                    pos[j] = np.array([atom.x, atom.y, atom.z])
                return pos
            
            # Add spring forces to atoms
            for i in range(1, len(NEB.states)-1):
                a = get_positions(NEB.states[i-1])
                b = get_positions(NEB.states[i])
                c = get_positions(NEB.states[i+1])
                
                real_force = np.array([np.empty([3]) for j in NEB.states[i]])
                for j,atom in enumerate(NEB.states[i]):
                    if j not in spring_atoms: continue
                    real_force[j] = np.array( [atom.fx,atom.fy,atom.fz] )

                # Find tangent
                tplus = c - b
                tminus = b - a
                dVmin = min( abs(V[i+1] - V[i]), abs(V[i-1]-V[i]) )
                dVmax = max( abs(V[i+1] - V[i]), abs(V[i-1]-V[i]) )

                if V[i+1] > V[i] and V[i] > V[i-1]:
                    tangent = tplus.copy()
                elif V[i+1] < V[i] and V[i] < V[i-1]:
                    tangent = tminus.copy()
                elif V[i+1] > V[i-1]:
                    tangent = tplus*dVmax + tminus*dVmin
                else:
                    tangent = tplus*dVmin + tminus*dVmax

                # Normalize tangent
                tangent_norm = np.sqrt(np.vdot(tangent,tangent))
                if tangent_norm != 0: tangent /= tangent_norm

                F_spring_parallel = k * (np.linalg.norm(tplus) - np.linalg.norm(tminus)) * tangent 
                
                F_real_perpendicular = real_force - (np.vdot(real_force,tangent)*tangent).reshape((-1,3))

                # Set NEB forces
                forces = F_spring_parallel + F_real_perpendicular

                for j,atom in enumerate(NEB.states[i]):
                    if j not in spring_atoms: continue
                    atom.fx, atom.fy, atom.fz = forces[j][0], forces[j][1], forces[j][2]

            # Remove net translation forces from the gradient
            #if fit_rigid is not None:
            if fit_rigid:
                net_translation_force = []
                for state in NEB.states[1:-1]:
                    net_force = np.zeros(3)
                    for a in state:
                        net_force += (a.fx, a.fy, a.fz)
                    net_translation_force.append(np.sqrt((net_force**2).sum()) / len(state))
                    for a in state:
                        a.fx -= net_force[0] / len(state)
                        a.fy -= net_force[1] / len(state)
                        a.fz -= net_force[2] / len(state)
                max_translation_force = units.convert("Ha/Ang","eV/Ang",max(net_translation_force))
            else:
                max_translation_force = 0
            
            # Set gradient
            NEB.gradient = []
            for state in NEB.states[1:-1]:
                for a in state:
                    NEB.gradient += [-a.fx, -a.fy, -a.fz] # Gradient of NEB.error

            # DEBUGGING CODE HERE
            #f_chk = []
            #for state in NEB.states[1:-1]:
                #for a in state:
                    #f_chk += [a.fx, a.fy, a.fz]
            #for i,f in enumerate(f_chk): f_chk[i] = units.convert_energy("Ha","eV",f)
            #f_chk = np.array(f_chk).reshape((-1,3))
            #print "Force Check 2", f_chk[3]
            #sys.exit()

            # Calculate RMS Force and Max force
            force_mags = [(a.fx**2+a.fy**2+a.fz**2)**0.5 for state in states[1:-1] for a in state]
            RMS_force = (sum([f**2 for f in force_mags])/float(len(force_mags)))**0.5
            NEB.RMS_force = RMS_force
            MAX_force = max(force_mags)
            NEB.MAX_force = MAX_force
                        
            # Print data
            V = V[:1] + [ units.convert_energy("Ha","kT_300",e-V[0]) for e in V[1:] ]
            MAX_energy = max(V)
            if NEB.prv_RMS == None or NEB.prv_RMS > RMS_force: rms = utils.color_set(float("%.4f" % units.convert_energy("Ha","eV",RMS_force)),'GREEN')
            else: rms = utils.color_set(float("%.4f" % units.convert_energy("Ha","eV",RMS_force)),'RED')
            if NEB.prv_MAX == None or NEB.prv_MAX > MAX_force: max_f = utils.color_set(float("%.4f" % units.convert_energy("Ha","eV",MAX_force)),'GREEN')
            else: max_f = utils.color_set(float("%.4f" % units.convert_energy("Ha","eV",MAX_force)),'RED')
            if NEB.prv_MAX_E == None or NEB.prv_MAX_E > MAX_energy: max_e = utils.color_set(float("%.1f" % MAX_energy),'GREEN')
            else: max_e = utils.color_set(float("%.1f" % MAX_energy),'RED')
            #print NEB.step, '%7.5g +' % V[0], ('%5.1f '*len(V[1:])) % tuple(V[1:]), rms, 'eV/Ang', max_f, 'eV/Ang', max_e, 'kT_300', ("%.4f" % max_translation_force), 'eV/Ang per atom'
            if NEB.step == 0:
                print("Step\tRMS_F (eV/Ang)\tMAX_F (eV/Ang)\tMAX_E (kT_300)\tMAX Translational Force (eV/Ang)\tEnergies (kT_300)\n----")
            print("%d\t%s\t\t%s\t\t%s\t\t%.4f" % (NEB.step, rms, max_f, max_e, max_translation_force)),
            print '    \t\t\t\t', '%7.5g +' % V[0], ('%5.1f '*len(V[1:])) % tuple(V[1:])
            
            sys.stdout.flush()
            
            if NEB.prv_RMS is None: NEB.prv_RMS = RMS_force
            NEB.prv_RMS = min(RMS_force, NEB.prv_RMS)

            if NEB.prv_MAX is None: NEB.prv_MAX = MAX_force
            NEB.prv_MAX = min(MAX_force, NEB.prv_MAX)

            if NEB.prv_MAX_E is None: NEB.prv_MAX_E = MAX_energy
            NEB.prv_MAX_E = min(MAX_energy, NEB.prv_MAX_E)

            # Set error
            #NEB.error = max(V)
            NEB.error = RMS_force
            #NEB.error = MAX_force
            #NEB.error = 0.2*max(V) + 0.8*RMS_force

            #NEB.error = 0.5*sum([a.fx**2+a.fy**2+a.fz**2 for state in states[1:-1] for a in state]) # Max(V) # Sum(energies)

            # Increment step
            NEB.step += 1
        
        # Calculate the "Error". This is what we minimize in our optimization algorithm.
        @staticmethod
        def get_error(coords):
            if NEB.error is None:
                NEB.calculate(coords)
            error = NEB.error
            NEB.error = None
            return error
        
        # Calculate the gradient.  This is negative of the force.
        @staticmethod
        def get_gradient(coords):
            if NEB.gradient is None:
                NEB.calculate(coords)
            gradient = NEB.gradient
            NEB.gradient = None #set to None so it will recalculate next time
            return np.array(gradient)

    def align_coordinates(r, B=None, H=None, return_matrix=False):
        # Prevent rotation or translation
        coord_count = 0
        st = NEB.states
        for s in st[1:-1]:
            for a in s:
                a.x, a.y, a.z = r[coord_count], r[coord_count+1], r[coord_count+2]
                coord_count += 3

        # Translate and rotate each frame to fit its neighbor
        # Note, procrustes will change st[-1] which is fine as we need this for spring
        # force calculations
        A = utils.procrustes(st)

        coord_count = 0
        for s in st[1:-1]:
            for a in s:
                r[coord_count:coord_count+3] = [a.x, a.y, a.z]
                coord_count += 3

        C = []
        R = block_diag(*A[0:-len(st[0])])
        if B is not None:
            for b in B:
                C.append(np.dot(b,R))
        if H is not None:
            # Note, to transform the Hessian matrix, it's not like a normal vector (as above)
            #H = R.T*H*R
            H = R*H*R.T

        if return_matrix:
            return A
        if B is None and H is None:
            return r
        elif B is not None and H is None:
            return r, C
        elif H is not None and B is None:
            return r, H
        else:
            return r, C, H

    def vproj(v1, v2):
        """
        Returns the projection of v1 onto v2
        Parameters:
            v1, v2: np vectors
        """
        mag2 = np.linalg.norm(v2)
        if mag2 == 0:
            print "Can't project onto a zero vector"
            return v1
        return v2 * (np.dot(v1, v2) / mag2)

    ######################################################################################
    ######################################################################################
     #######  ########  ######## #### ##     ## #### ######## ######## ########   ######  
    ##     ## ##     ##    ##     ##  ###   ###  ##       ##  ##       ##     ## ##    ## 
    ##     ## ##     ##    ##     ##  #### ####  ##      ##   ##       ##     ## ##       
    ##     ## ########     ##     ##  ## ### ##  ##     ##    ######   ########   ######  
    ##     ## ##           ##     ##  ##     ##  ##    ##     ##       ##   ##         ## 
    ##     ## ##           ##     ##  ##     ##  ##   ##      ##       ##    ##  ##    ## 
     #######  ##           ##    #### ##     ## #### ######## ######## ##     ##  ######  
    ######################################################################################
    ######################################################################################

    def steepest_descent(f, r, fprime, alpha=0.05, maxiter=1000, gtol=1E-3, fit_rigid=True): #better, but tends to push error up eventually, especially towards endpoints.
        step = 0
        while (NEB.RMS_force > gtol) and (step < maxiter):
            if NEB.FMAX is not None and NEB.MAX_force < NEB.FMAX: return FMAX_CONVERGENCE

            f = np.array(-fprime(r)).reshape((-1,3))

            max_step_length = np.sqrt(((f)**2).sum(axis=1).max())
            # Scale if largest step to be taken will become larger than alpha:
            if max_step_length > 1.0:
                dr = f * alpha / max_step_length
            else:
                dr = f * alpha

            r += dr.flatten()
            if fit_rigid: r = align_coordinates(r)
            step += 1

        if NEB.RMS_force <= gtol: return GTOL_CONVERGENCE
        if step >= maxiter: return MAXITER_CONVERGENCE
        return FAIL_CONVERGENCE

    def quick_min_optimizer(f, r, nframes, fprime, dt=0.1, step_max=0.1, euler=False, viscosity=0.1, maxiter=1000, gtol=1E-3, fit_rigid=True): # dt = fs, step_max = angstroms, viscosity = 1/fs
        v = np.array([0.0 for x in r])
        acc = np.array([0.0 for x in r])

        masses = []
        for s in states[1:-1]:
            for a in s:
                m = units.elem_weight(a.element)
                masses += [m, m, m]
        masses = np.array(masses)

        step = 0
        while (NEB.RMS_force > gtol) and (step < maxiter):
            if NEB.FMAX is not None and NEB.MAX_force < NEB.FMAX: return FMAX_CONVERGENCE
            forces = -fprime(r) # Get the forces

            # Get the parallel velocity if it's in the same direction as the force
            # Note, this must be frame based, not individual atom based or total based
            natoms = len(v)/(3*(nframes-2))
            #print("\nNum Atoms = %lg\n" % natoms)
            for i in range(1,nframes-1):
                #print("Frame %d of %d:" % (i,nframes)),
                low = (i-1)*natoms*3
                high = i*natoms*3

                force_direction = forces[low:high]/np.linalg.norm(forces[low:high]) # Get the direction of the force
                
                import random
                
                if np.dot(forces[low:high],v[low:high]) > 0.0:
                    #print sum(v[low:high]),
                    v[low:high] = vproj(v[low:high],forces[low:high])
                    #print sum(v[low:high])
                else:
                    v[low:high] *= 0.0
                    print 'zeroed velocities in frame %d' % i
                
                speed = np.linalg.norm(v[low:high])
                if speed*dt > step_max:
                    max_speed = step_max/dt
                    v[low:high] *= max_speed / speed
                
            if euler:
                #make Euler step
                v += dt * forces

                #limit distance moved
                #for i in range(len(v)):
                #   if v[i]*dt > step_max: v[i] = step_max/dt
                #   if v[i]*dt <-step_max: v[i] =-step_max/dt

                for i in range(1,nframes-1):
                    low = (i-1)*natoms*3
                    high = i*natoms*3
                    speed = np.linalg.norm(v[low:high])
                    if speed*dt > step_max:
                        max_speed = step_max/dt
                        v[low:high] *= max_speed / speed

                #move atoms
                r += v * dt
            else:
                #make Verlet step
                a_new = forces/masses
                a_new -= v*viscosity
                
                dx = v*dt + 0.5*acc*dt**2
                #limit distance moved
                #for i in range(len(r)):
                #   if dx[i] > step_max: dx[i] = step_max
                #   if dx[i] <-step_max: dx[i] =-step_max
                
                r_new = r + dx
                v_new = v + (acc + a_new)*0.5 * dt
                r = r_new
                v = v_new
                acc = a_new
            
            if fit_rigid: r = align_coordinates(r)

            step += 1
        if NEB.RMS_force <= gtol: return GTOL_CONVERGENCE
        if step >= maxiter: return MAXITER_CONVERGENCE
        return FAIL_CONVERGENCE

    def fire_optimizer(f, r, nframes, fprime, dt = 0.1, dtmax = 1.0, step_max = 0.2, maxiter=1000, gtol=1E-3,
                        Nmin = 5, finc = 1.1, fdec = 0.5, astart = 0.1, fa = 0.99, euler = True, fit_rigid=True):

        v = np.array([0.0 for x in r])
        Nsteps = 0
        acc = astart

        step = 0
        while (NEB.RMS_force > gtol) and (step < maxiter):
            if NEB.FMAX is not None and NEB.MAX_force < NEB.FMAX: return FMAX_CONVERGENCE
            # Get forces and number of atoms
            forces = -fprime(r)
            natoms = len(v)/(3*(nframes-2))
            
            if np.dot(v,forces) > 0.0:
                # If velocity in direction of forces, speed up
                v = (1.0-acc)*v + acc*np.linalg.norm(v)*(forces/np.linalg.norm(forces))
                if(Nsteps>Nmin):
                    dt = min(dt*finc,dtmax)
                    acc *= fa
                Nsteps += 1
            else:
                # If not, slow down
                v *= 0.0
                acc = astart
                dt *= fdec
                Nsteps = 0

            if euler:
                #make Euler step
                v += dt * forces

                #limit velocities
                for i in range(len(v)):
                    if v[i]*dt > step_max: v[i] = step_max/dt
                    if v[i]*dt <-step_max: v[i] =-step_max/dt

                #move atoms
                r += v * dt

            if fit_rigid: r = align_coordinates(r)

            step += 1

        if NEB.RMS_force <= gtol: return GTOL_CONVERGENCE
        if step >= maxiter: return MAXITER_CONVERGENCE
        return FAIL_CONVERGENCE

    def lbfgs_optimizer(target_function, initial_coordinates, target_gradient,
            step_size=0.1, step_size_adjustment=0.5, armijio_line_search_factor=1E-4, reset_when_in_trouble=True, linesearch='armijo',
            gradient_tolerance=1E-3, max_iterations=1000, reset_step_size=reset,
            MAX_STEP=0.2, fit_rigid=True,
            display=0, callback=None, max_steps_remembered=Nmax):

        # These are values deemed good for DFT NEB and removed from parameter space for simplification
        MIN_STEP=1E-8
        BACKTRACK_EPS=1E-3

        if display > 2:
            print("\nValues in bfgs_optimize code:")
            print("\tstep size = %lg, step size adjustment = %lg, reset_when_in_trouble = %s" % (step_size, step_size_adjustment, str(reset_when_in_trouble)))
            print("\tgtol = %lg, max_iterations = %d, MAX_STEP = %lg" % (gradient_tolerance, max_iterations, MAX_STEP))
            print("\t-------------------")
            print("\treset_step_size = %s, MIN_STEP = %lg, BACKTRACK_EPS = %lg" % (str(reset_step_size), MIN_STEP, BACKTRACK_EPS))
            print("\tfit_rigid = %s\n" % str(fit_rigid))
       
        function_call_counter = 0

        # Set max_iterations if not set
        if max_iterations is None:
            max_iterations = 200 * len(initial_coordinates)

        # Ensure coordinates are in the correct format
        initial_coordinates = np.asarray(initial_coordinates).flatten()
        if initial_coordinates.ndim == 0:
            initial_coordinates.shape = (1,)
        if fit_rigid: current_coordinates = align_coordinates(initial_coordinates)
        else: current_coordinates = copy.deepcopy(initial_coordinates)

        # Initialize stored coordinates and gradients
        stored_coordinates = []
        stored_gradients = []

        # Get gradient and store your old func_max
        current_gradient = target_gradient(current_coordinates)
        if target_function is not None:
            old_fval = target_function(current_coordinates)
        function_call_counter += 1

        # Hold original values
        ALPHA_CONST = step_size
        BETA_CONST = step_size_adjustment
        RESET_CONST = reset_step_size

        # Get function to describe linesearch
        if linesearch is 'armijo':
            if display > 1: print("armijo linesearch "),
            def check_backtrack(f1,f0,gk,pk,armijio_line_search_factor,step_size):
                return f1-f0 > armijio_line_search_factor*step_size*np.dot(gk,pk)
        else:
            if display > 1: print("default linesearch "),
            def check_backtrack(f1,f0,pk,gk,armijio_line_search_factor,step_size):
                return (f1-f0)/(abs(f1)+abs(f0)) > BACKTRACK_EPS

        backtrack, loop_counter, warnflag = 0, 0, 0
        while (NEB.RMS_force > gradient_tolerance) and (function_call_counter < max_iterations):
            if NEB.FMAX is not None and NEB.MAX_force < NEB.FMAX: return FMAX_CONVERGENCE

            if display > 1:
                print("Step %d, " % loop_counter),

            # Get your step direction
            #step_direction = -np.dot(current_Hessian, current_gradient)
            
            def BFGS_multiply(s,y,grad): #http://aria42.com/blog/2014/12/understanding-lbfgs/, with corrections
                r = copy.deepcopy(grad)
                indices = xrange(len(s))
                #compute right product
                alpha = np.zeros( len(s) )
                for i in reversed(indices):
                    rho_i = 1.0 / np.dot(y[i],s[i])
                    alpha[i] = rho_i * np.dot(s[i],r)
                    r -= alpha[i] * y[i]
                
                # Only multiply by approximate inv_hessian if we have stored coordinates
                if len(s) > 0:
                    r *= np.dot(y[-1],s[-1]) / np.dot(y[-1],y[-1])

                #compute left product
                for i in indices:
                    rho_i = 1.0 / np.dot(y[i],s[i])
                    beta = rho_i * np.dot(y[i],r)
                    r += ( alpha[i] - beta )*s[i]
                return r
            
            step_direction = -BFGS_multiply(list(reversed(stored_coordinates)), list(reversed(stored_gradients)), current_gradient)
            
            # Renorm to remove the effect of Hessian not being unit
            i = 0
            while i < len(step_direction):
                # Get the distance the atom will move
                a,b,c = step_direction[i],step_direction[i+1],step_direction[i+2]
                chk = float((a**2+b**2+c**2)**0.5)
                a,b,c = current_gradient[i],current_gradient[i+1],current_gradient[i+2]
                scale = float((a**2+b**2+c**2)**0.5)
                
                step_direction[i] *= (scale / chk)
                step_direction[i+1] *= (scale / chk)
                step_direction[i+2] *= (scale / chk)
                i += 3

            # If we are doing unreasonably small step sizes, quit
            if abs(max(step_direction*step_size)) < MIN_STEP:
                if display > 1:
                    print("Error - Step size unreasonable (%lg)" 
                                % abs(max(step_direction*step_size))),
                warnflag = 2
                break

            # If we have too large of a step size, set to max
            # Loop through atoms
            i, max_step_flag = 0, False
            while i < len(step_direction):
                # Get the distance the atom will move
                a,b,c = step_direction[i]*step_size,step_direction[i+1]*step_size,step_direction[i+2]*step_size
                chk = float((a**2+b**2+c**2)**0.5)

                # If d = sqrt(a^2+b^2+c^2) > MAX_STEP, scale by MAX_STEP/d
                if chk > MAX_STEP:
                    max_step_flag = True
                    step_direction[i] *= (MAX_STEP / chk)
                    step_direction[i+1] *= (MAX_STEP / chk)
                    step_direction[i+2] *= (MAX_STEP / chk)
                i += 3

            # As we are changing values manually, this is no longer
            # the BFGS(Hess) algorithm so reset the Inverse Hessian
            if max_step_flag:
                if display > 1:
                    print("Warning - Setting step to max step size"),
                if reset_when_in_trouble:
                    stored_coordinates = []
                    stored_gradients = []

            # Hold new parameters
            new_coordinates = current_coordinates + step_size * step_direction

            #if fit_rigid:
            #    rotation = align_coordinates(new_coordinates, return_matrix=True)
            #    print rotation
            #    for i in xrange(len(stored_coordinates)):
            #        for j in range(0, len(new_coordinates), 3):
            #            stored_coordinates[i][j], stored_coordinates[i][j+1], stored_coordinates[i][j+2] = np.dot((stored_coordinates[i][j], stored_coordinates[i][j+1], stored_coordinates[i][j+2]), rotation)
            #            stored_gradients[i][j], stored_gradients[i][j+1], stored_gradients[i][j+2] = np.dot((stored_gradients[i][j], stored_gradients[i][j+1], stored_gradients[i][j+2]), rotation)
            
            # Get the new gradient
            new_gradient = target_gradient(new_coordinates)

            # Check if max has increased
            if target_function is not None:
                fval = target_function(new_coordinates)
            function_call_counter += 1

            if target_function is not None and check_backtrack(fval, old_fval, new_gradient, step_direction, armijio_line_search_factor, step_size):
                # Step taken overstepped the minimum.  Lowering step size
                if display > 1:
                    print("\tResetting System as %lg > %lg!"
                            % (fval, old_fval))
                    print("\talpha: %lg" % step_size),

                step_size *= np.float64(step_size_adjustment)

                if display > 1:
                    print("-> %lg\n" % step_size)

                # Reset the Inverse Hessian if desired.
                # It is still up for debate if this is to be recommended or not.  As the 
                # inverse hessian corrects itself, it might not be important to do this.
                if reset_when_in_trouble:
                    stored_coordinates = []
                    stored_gradients = []
                backtrack += 1
                reset_step_size = RESET_CONST
                continue

            # This allows for the edge case in which after decreasing step_size, a situation arises
            # in which larger alphas are acceptable again. Thus, we reset to the original step_size
            elif reset_step_size is not None:
                reset_step_size -= 1
                # If we want to reset_step_size and step_size has been decreased before, set to initial vals
                if reset_step_size < 0 and step_size < ALPHA_CONST:
                    if display > 1:
                        print("\tResetting Alpha, Beta, Reset and Inverse Hessian")
                    step_size, step_size_adjustment, reset_step_size = ALPHA_CONST, BETA_CONST, RESET_CONST
                    # Once again, debatable if we want this here.  When reseting step sizes we
                    # might have a better H inverse than the Identity would be.
                    if reset_when_in_trouble:
                        stored_coordinates = []
                        stored_gradients = []
                    continue
                # If we want to reset_step_size and we've never decreased before, we can take larger steps
                # We increase step sizes by a factor of 1/step_size_adjustment
                elif reset_step_size < 0 and step_size >= ALPHA_CONST:
                    if display > 1:
                        print("\tIncreasing step size: %lg ->" % step_size),
                    step_size /= step_size_adjustment
                    if display > 1:
                        print("%lg,\t" % step_size),
            
            # Recalculate change_in_coordinates to maintain the secant condition
            change_in_coordinates = new_coordinates - current_coordinates
            
            # Store new max value in old_max for future comparison
            if target_function is not None:
                old_fval = fval

            # Get difference in gradients for further calculations
            change_in_gradient = new_gradient - current_gradient

            #store past results to build up curvature information
            stored_coordinates.append( change_in_coordinates )
            stored_gradients.append( change_in_gradient )
            
            if len(stored_coordinates)>max_steps_remembered:
                stored_coordinates = stored_coordinates[1:]
                stored_gradients = stored_gradients[1:]

            try:  # this was handled in numeric, let it remain for more safety
                rhok = 1.0 / (np.dot(change_in_gradient, change_in_coordinates))
            except ZeroDivisionError:
                rhok = 1000.0
                if display > 1:
                    print("Divide-by-zero encountered: rhok assumed large")
            if np.isinf(rhok):  # this is patch for np
                rhok = 1000.0
                if display > 1:
                    print("Divide-by-zero encountered: rhok assumed large")


            # Run BFGS Update for the Inverse Hessian
            #A1 = I - change_in_coordinates[:, np.newaxis] * change_in_gradient[np.newaxis, :] * rhok
            #A2 = I - change_in_gradient[:, np.newaxis] * change_in_coordinates[np.newaxis, :] * rhok
            #current_Hessian = np.dot(A1, np.dot(current_Hessian, A2)) + (rhok * change_in_coordinates[:, np.newaxis] * change_in_coordinates[np.newaxis, :])

            if display > 1:
                print("fval %lg" % (fval))

            # Store new parameters, as it has passed the check
            current_coordinates = new_coordinates
            current_gradient = new_gradient

            # If callback is desired
            if callback is not None:
                callback(current_coordinates)

            # Increment the loop counter
            loop_counter += 1

            if (NEB.RMS_force <= gradient_tolerance):
                break
            
        if target_function is not None:
            fval = old_fval
        else:
            fval = float('NaN')

        if np.isnan(fval):
            # This can happen if the first call to f returned NaN;
            # the loop is then never entered.
            warnflag = 2

        if warnflag == 2:
            if display == 1:
                print("Warning: Loss of precision.")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % function_call_counter)

        elif loop_counter >= max_iterations:
            warnflag = 1
            if display == 1:
                print("Warning: Maximum Iteration was exceeded.")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % function_call_counter)
        else:
            if display == 1:
                print("Success!")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % function_call_counter)

        if NEB.RMS_force <= gtol: return GTOL_CONVERGENCE
        if loop_counter >= maxiter: return MAXITER_CONVERGENCE
        return FAIL_CONVERGENCE


######################################################################################################
######################################################################################################
######################################################################################################
# ORIGINAL BFGS CODE
######################################################################################################
######################################################################################################
######################################################################################################

    def bfgs_optimizer(target_function, initial_coordinates, target_gradient,
            step_size=0.1, step_size_adjustment=0.5, armijio_line_search_factor=1E-4, reset_when_in_trouble=True, linesearch='armijo',
            gradient_tolerance=1E-3, max_iterations=1000, reset_step_size=reset,
            MAX_STEP=0.2, fit_rigid=True,
            display=0, callback=None):

        # Wrap functions together
        def target(coordinates, function_call_counter):
            function_call_counter += 1
            if target_function is not None:
                return target_function(coordinates), target_gradient(coordinates), function_call_counter
            else:
                return None, target_gradient(coordinates), function_call_counter

        # These are values deemed good for DFT NEB and removed from parameter space for simplification
        MIN_STEP=1E-8
        BACKTRACK_EPS=1E-3
        function_call_counter = 0        

        if display > 2:
            print("\nValues in bfgs_optimize code:")
            print("\tstep size = %lg, step size adjustment = %lg, reset_when_in_trouble = %s" % (step_size, step_size_adjustment, str(reset_when_in_trouble)))
            print("\tgtol = %lg, max_iterations = %d, MAX_STEP = %lg" % (gradient_tolerance, max_iterations, MAX_STEP))
            print("\t-------------------")
            print("\treset_step_size = %s, MIN_STEP = %lg, BACKTRACK_EPS = %lg" % (str(reset_step_size), MIN_STEP, BACKTRACK_EPS))
            print("\tfit_rigid = %s\n" % str(fit_rigid))

        # Ensure coordinates are in the correct format
        initial_coordinates = np.asarray(initial_coordinates).flatten()
        if initial_coordinates.ndim == 0:
            initial_coordinates.shape = (1,)
        current_coordinates = copy.deepcopy(initial_coordinates)

        # Initialize inv Hess and Identity matrix
        I = np.eye(len(current_coordinates), dtype=int)
        current_Hessian = I

        # Get gradient and store your old func_max
        old_fval, current_gradient, function_call_counter = target(current_coordinates, function_call_counter)

        # Hold original values
        ALPHA_CONST = step_size
        BETA_CONST = step_size_adjustment
        RESET_CONST = reset_step_size

        # Get function to describe linesearch
        if linesearch is 'armijo':
            if display > 1: print("armijo linesearch "),
            def check_backtrack(f1,f0,gk,pk,armijio_line_search_factor,step_size):
                return f1-f0 > armijio_line_search_factor*step_size*np.dot(gk,pk)
        else:
            if display > 1: print("default linesearch "),
            def check_backtrack(f1,f0,pk,gk,armijio_line_search_factor,step_size):
                return (f1-f0)/(abs(f1)+abs(f0)) > BACKTRACK_EPS

        backtrack, loop_counter, warnflag = 0, 0, 0
        while (NEB.RMS_force > gradient_tolerance) and (function_call_counter < max_iterations):
            if NEB.FMAX is not None and NEB.MAX_force < NEB.FMAX: return FMAX_CONVERGENCE

            if display > 1:
                print("Trace of current_Hessian = %lg" % float(np.matrix.trace(current_Hessian)))
                print("Step %d, " % loop_counter),

            # Get your step direction and renorm to remove the effect of current_Hessian not being unit
            step_direction = -np.dot(current_Hessian, current_gradient).reshape((-1,3))
            force_mags = (current_gradient.reshape((-1,3))**2).sum(axis=1)
            scalar = np.sqrt(force_mags / (step_direction**2).sum(axis=1))
            step_direction = (step_direction.T * scalar).T

            # If we are doing unreasonably small step sizes, quit
            step_lengths = np.sqrt( (step_direction**2).sum(axis=1) ) * step_size
            if max(step_lengths) < MIN_STEP:
                if display > 1: print("Error - Step size unreasonable (%lg)" % abs(max(step_length))),
                warnflag = 2
                break

            # If we have too large of a step size, set to max
            indices_of_large_steps = [(i,s) for i,s in enumerate(step_lengths) if s > MAX_STEP]
            for i,s in indices_of_large_steps: step_direction[i] *= MAX_STEP/s
            max_step_flag = len(indices_of_large_steps) > 0
            step_direction = step_direction.flatten()

            # As we are changing values manually, this is no longer
            # the BFGS(Hess) algorithm so reset the Inverse Hessian
            # -> This is because we're no longer on the eigendirection
            if max_step_flag:
                if display > 1: print("Warning - Setting step to max step size"),
                if reset_when_in_trouble: current_Hessian = I.copy()

            # Hold new parameters
            # We do so because we may not want to keep them if it makes the next step bad.
            new_coordinates = current_coordinates + step_size * step_direction
            if fit_rigid:
                new_coordinates, C, current_Hessian_tmp = align_coordinates(new_coordinates, [current_gradient, current_coordinates], current_Hessian)
                current_gradient_tmp, current_coordinates_tmp = C

            # Get the new gradient and check if max has increased
            fval, new_gradient, function_call_counter = target(new_coordinates, function_call_counter)

            if target_function is not None and check_backtrack(fval, old_fval, new_gradient, step_direction, armijio_line_search_factor, step_size):
                # Step taken overstepped the minimum.  Lowering step size
                if display > 1:
                    print("\tResetting System as %lg > %lg!" % (fval, old_fval))
                    print("\talpha: %lg" % step_size),

                step_size *= np.float64(step_size_adjustment)

                if display > 1: print("-> %lg\n" % step_size)

                # Reset the Inverse Hessian if desired.
                # It is still up for debate if this is to be recommended or not.  As the 
                # inverse hessian corects itself, it might not be important to do this.
                if reset_when_in_trouble: current_Hessian = I
                backtrack += 1
                reset_step_size = RESET_CONST
                continue

            # This allows for the edge case in which after decreasing step_size, a situation arises
            # in which larger alphas are acceptable again. Thus, we reset to the original step_size
            elif reset_step_size is not None:
                reset_step_size -= 1
                # If we want to reset_step_size and step_size has been decreased before, set to initial vals
                if reset_step_size < 0 and step_size < ALPHA_CONST:
                    if display > 1: print("\tResetting Alpha, Beta, Reset and Inverse Hessian")
                    step_size, step_size_adjustment, reset_step_size = ALPHA_CONST, BETA_CONST, RESET_CONST
                    # Once again, debatable if we want this here.  When reseting step sizes we
                    # might have a better H inverse than the Identity would be.
                    if reset_when_in_trouble: current_Hessian = I
                    continue
                # If we want to reset_step_size and we've never decreased before, we can take larger steps
                # We increase step sizes by a factor of 1/step_size_adjustment
                elif reset_step_size < 0 and step_size >= ALPHA_CONST:
                    if display > 1:
                        print("\tIncreasing step size: %lg ->" % step_size),
                    step_size /= step_size_adjustment
                    if display > 1:
                        print("%lg,\t" % step_size),
            
            # If the step was good, we want to store the rotated values
            if fit_rigid:
                current_gradient, current_coordinates, current_Hessian = current_gradient_tmp, current_coordinates_tmp, current_Hessian_tmp

            # Recalculate change_in_coordinates to maintain the secant condition
            change_in_coordinates = new_coordinates - current_coordinates
            
            # Store new max value in old_max for future comparison
            if target_function is not None: old_fval = fval

            # Get difference in gradients for further calculations
            change_in_gradient = new_gradient - current_gradient

            try:  # this was handled in numeric, let it remaines for more safety
                rhok = 1.0 / (np.dot(change_in_gradient, change_in_coordinates))
            except ZeroDivisionError:
                rhok = 1000.0
                if display > 1:
                    print("Divide-by-zero encountered: rhok assumed large")
            if np.isinf(rhok):  # this is patch for np
                rhok = 1000.0
                if display > 1:
                    print("Divide-by-zero encountered: rhok assumed large")


            # Run BFGS Update for the Inverse Hessian
            A1 = I - change_in_coordinates[:, np.newaxis] * change_in_gradient[np.newaxis, :] * rhok
            A2 = I - change_in_gradient[:, np.newaxis] * change_in_coordinates[np.newaxis, :] * rhok
            current_Hessian = np.dot(A1, np.dot(current_Hessian, A2)) + (rhok * change_in_coordinates[:, np.newaxis] * change_in_coordinates[np.newaxis, :])

            if display > 1: print("fval %lg" % (fval))

            # Store new parameters, as it has passed the check
            current_coordinates = new_coordinates
            current_gradient = new_gradient

            # If callback is desired
            if callback is not None:
                callback(current_coordinates)

            # Increment the loop counter
            loop_counter += 1

            if (NEB.RMS_force <= gradient_tolerance):
                break
            
        if target_function is not None: fval = old_fval
        else: fval = float('NaN')

        # This can happen if the first call to f returned NaN;
        # the loop is then never entered.
        if np.isnan(fval): warnflag = 2

        if warnflag == 2:
            if display == 1:
                print("Warning: Loss of precision.")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % function_call_counter)

        elif loop_counter >= max_iterations:
            warnflag = 1
            if display == 1:
                print("Warning: Maximum Iteration was exceeded.")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % function_call_counter)
        else:
            if display == 1:
                print("Success!")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % function_call_counter)

        if NEB.RMS_force <= gtol: return GTOL_CONVERGENCE
        if loop_counter >= maxiter: return MAXITER_CONVERGENCE
        return FAIL_CONVERGENCE

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################



    #######################################################################################################################
    #######################################################################################################################
    ########  ##     ## ##    ##     ######  #### ##     ## ##     ## ##          ###    ######## ####  #######  ##    ## 
    ##     ## ##     ## ###   ##    ##    ##  ##  ###   ### ##     ## ##         ## ##      ##     ##  ##     ## ###   ## 
    ##     ## ##     ## ####  ##    ##        ##  #### #### ##     ## ##        ##   ##     ##     ##  ##     ## ####  ## 
    ########  ##     ## ## ## ##     ######   ##  ## ### ## ##     ## ##       ##     ##    ##     ##  ##     ## ## ## ## 
    ##   ##   ##     ## ##  ####          ##  ##  ##     ## ##     ## ##       #########    ##     ##  ##     ## ##  #### 
    ##    ##  ##     ## ##   ###    ##    ##  ##  ##     ## ##     ## ##       ##     ##    ##     ##  ##     ## ##   ### 
    ##     ##  #######  ##    ##     ######  #### ##     ##  #######  ######## ##     ##    ##    ####  #######  ##    ##
    #######################################################################################################################
    #######################################################################################################################

    n = NEB(name, states, theory, extra_section, k, fmax, start_from=start_from)

    # Output for user
    if opt == 'BROYDEN_ROOT':
        scipy.optimize.broyden1(NEB.get_gradient, np.array(NEB.coords_start), alpha=float(alpha), verbose=(disp != 0))
    elif opt == 'QM':
        conv_flag = quick_min_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.nframes, 
            fprime=NEB.get_gradient, dt=dt, viscosity=viscosity, step_max=step_max, euler=euler, maxiter=maxiter, gtol=gtol, fit_rigid=fit_rigid)
    elif opt == 'FIRE':
        conv_flag = fire_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.nframes, 
            fprime=NEB.get_gradient, dt=dt, dtmax=dtmax, step_max=step_max,
            Nmin=Nmin, finc=finc, fdec=fdec, astart=astart, fa=fa, euler=euler, maxiter=maxiter, gtol=gtol, fit_rigid=fit_rigid)
    elif opt == 'BFGS':
        conv_flag = bfgs_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.get_gradient,
            step_size=float(alpha), step_size_adjustment=float(beta), armijio_line_search_factor=float(tau), reset_when_in_trouble=H_reset,
            gradient_tolerance=float(gtol), max_iterations=int(maxiter), fit_rigid=fit_rigid, linesearch=linesearch, reset_step_size=reset,
            MAX_STEP=float(step_max), display=disp
            )
    elif opt == 'LBFGS':
        conv_flag = lbfgs_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.get_gradient,
            step_size=float(alpha), step_size_adjustment=float(beta), armijio_line_search_factor=float(tau), reset_when_in_trouble=H_reset,
            gradient_tolerance=float(gtol), max_iterations=int(maxiter), fit_rigid=fit_rigid, linesearch=linesearch, reset_step_size=reset,
            MAX_STEP=float(step_max), max_steps_remembered=Nmax, display=disp
            )
    elif opt == 'SD':
        conv_flag = steepest_descent(NEB.get_error, np.array(NEB.coords_start), fprime=NEB.get_gradient, alpha=alpha, maxiter=maxiter, gtol=gtol, fit_rigid=fit_rigid)
    else:
        print("\nERROR - %s optimizations method does not exist! Choose from the following:" % str(opt))
        print("\t1. BFGS")
        print("\t2. LBFGS")
        print("\t3. QM")
        print("\t4. SD")
        print("\t5. FIRE")
        print("\t6. BROYDEN_ROOT\n")
        sys.exit()

    if conv_flag == FAIL_CONVERGENCE:
        print("\nNEB failed to converge.")
    elif conv_flag == MAXITER_CONVERGENCE:
        print("\nNEB quit after reaching the specified maximum number of iterations.")
    elif conv_flag == FMAX_CONVERGENCE:
        print("\nNEB converged the maximum force.")
    elif conv_flag == GTOL_CONVERGENCE:
        print("\nNEB converged the RMS force.")