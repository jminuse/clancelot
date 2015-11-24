from merlin import *
import re, shutil, copy

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

def neb(name, states, theory, extra_section='', spring_atoms=None, procs=1, queue=None,
        disp=0, fit_rigid=True, center=True, k=0.1837,
        DFT='g09', opt='BFGS', gtol=1e-5, maxiter=1000,
        alpha=0.1, beta=0.6, tau=1E-3, reset=60, H_reset=True,
        viscosity=0.1, dtmax=1.0, Nmin=5, finc=1.1, fdec=0.5, astart=0.1, fa=0.99,
        step_min=1E-10, step_max=0.2, bt_max=None, linesearch='backtrack', L2norm=True, bt_eps=1E-3,
        dt = 0.1, euler=True, force=True, mem=25, blurb=None, initial_guess=None): 
    
    # If using test code, import path so we import correct scipy.optimize.
    if opt=='BFGS2': sys.path.insert(1,'/fs/home/hch54/scipy_mod/scipy/')
    import scipy.optimize
    import numpy as np

    DFT = DFT.lower().strip()

    # Contemplating a force addition of NoSymm to route if (1) DFT='g09' and (2) NoSymm not said

    # Set which atoms will be affected by virtual springs
    if not spring_atoms: # If not given, select all
        spring_atoms = range(len(states[0]))
    elif type(spring_atoms)==str: # A list of element names
        elements = spring_atoms.split()
        spring_atoms = [i for i,a in enumerate(states[0]) if a.element in elements]

    # Output for user
    if opt == 'BROYDEN_ROOT':
        print("\nRunning neb with optimizaiton method %s" % str(opt))   
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
    elif opt in ['BFGS','BFGS2']:
        print("\nRunning neb with optimization method %s" % str(opt))
        print("\talpha = %lg, beta = %lg" % (alpha, beta)),
        if linesearch == 'armijo': print(", tau = %lg" % tau)
        else: print("")
        print("\tH_reset = %s, reset = %s, linesearch = %s" % (str(H_reset), str(reset), linesearch))
        print("\tstep_min = %lg, step_max = %lg, L2norm = %s" % (step_min, step_max, str(L2norm)))
        print("\tbt_max = %s, bt_eps = %lg" % (str(bt_max), bt_eps))
    elif opt == 'SD':
        print("\nRunning neb with optimization method %s" % str(opt))
        print("\talpha = %lg" % alpha)
    else:
        print("\nERROR - %s optimizations method does not exist! Choose from the following:" % str(opt))
        print("\t1. BFGS")
        print("\t2. BFGS2")
        print("\t3. QM")
        print("\t4. SD")
        print("\t5. FIRE")
        print("\t6. BROYDEN_ROOT\n")
        sys.exit()

    print("Spring Constant for NEB: %lg Ha/Ang" % k)
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
        def __init__(self, name, states, theory, extra_section='', k=0.1837, fit_rigid=True):
            NEB.name = name
            NEB.states = states
            NEB.theory = theory
            NEB.extra_section = extra_section
            NEB.k = k
            NEB.prv_RMS = None
            NEB.convergence_criteria = gtol
            NEB.convergence = float('inf')
            NEB.nframes = len(states)

            if fit_rigid: 
                utils.procrustes(NEB.states) # Fit rigid before relaxing
            
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
                
                # Set the extra section to be whatever specified originally
                extra_section = NEB.extra_section.strip()
                # If g09, ensure two newlines at end
                if DFT=='g09':
                    extra_section += '\n\n'
                
                if NEB.step>0:
                    if (i==0 or i==len(NEB.states)-1): # Only use DFT on endpoints on first step, because they don't change
                        continue
                    if DFT=='g09':
                        guess = ' Guess=Read'
                    elif DFT=='orca':
                        guess = ' MOREAD'
                        tmp = '%%moinp "../%s-%d-%d/%s-%d-%d.orca.gbw"' % (NEB.name,NEB.step-1,i,NEB.name,NEB.step-1,i)
                        extra_section = tmp + extra_section.strip()
                else:
                    if initial_guess:
                        if DFT=='g09':
                            guess = ' Guess=Read'
                        elif DFT=='orca':
                            guess = ' MOREAD'
                            tmp = '%%moinp "../%s/%s.orca.gbw\n"' % (initial_guess,initial_guess)
                            extra_section = tmp + extra_section.strip()
                    else:
                        guess = '' #no previous guess for first step

                if DFT=='g09':
                    running_jobs.append( g09.job('%s-%d-%d'%(NEB.name,NEB.step,i), NEB.theory+' Force'+guess, state, procs=procs, queue=queue, force=force, previous=('%s-%d-%d'%(NEB.name,NEB.step-1,i)) if NEB.step>0 else initial_guess, extra_section=extra_section, neb=[True,'%s-%%d-%%d'%(NEB.name),len(NEB.states),i], mem=mem) )
                elif DFT=='orca':
                	running_jobs.append( orca.job('%s-%d-%d'%(NEB.name,NEB.step,i), NEB.theory+guess, state, extra_section=extra_section, grad=True, procs=procs, queue=queue) )

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

                if DFT=='g09':
                    result = g09.parse_atoms('%s-%d-%d' % (NEB.name, step_to_use, i), check_convergence=False, parse_all=False)
                    if not result:
                        raise Exception('parse_atoms failed')
                    new_energy, new_atoms = result
                elif DFT=='orca':
                    new_atoms, new_energy = orca.engrad_read('%s-%d-%d' % (NEB.name, step_to_use, i))

                energies.append(new_energy)

                if DFT == 'g09':
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
                            a.fx = b.fx; a.fy = b.fy; a.fz = b.fz
                elif DFT == 'orca':
                    for a,b in zip(state, new_atoms):
						a.fx = b.fx; a.fy = b.fy; a.fz = b.fz
            # V = potential energy from DFT. energies = V+springs
            V = copy.deepcopy(energies) 
            # Reset convergence check
            NEB.convergence = 0

            # Add spring forces to atoms
            for i,state in enumerate(NEB.states):
                if i==0 or i==len(NEB.states)-1: continue # Don't change first or last state
                for j,b in enumerate(state):
                    if j in spring_atoms:
                        a,c = NEB.states[i-1][j], NEB.states[i+1][j]
                    
                        # Find tangent
                        tplus = np.array( [ c.x-b.x, c.y-b.y, c.z-b.z ] )
                        tminus = np.array( [ b.x-a.x, b.y-a.y, b.z-a.z ] )
                        dVmin = min(abs(V[i+1]-V[i]), abs(V[i-1]-V[i]))
                        dVmax = max(abs(V[i+1]-V[i]), abs(V[i-1]-V[i]))
                        if V[i+1]>V[i] and V[i]>V[i-1]: # Not at an extremum, trend of V is up
                            tangent = tplus
                        elif V[i+1]<V[i] and V[i]<V[i-1]: # Not at an extremum, trend of V is down
                            tangent = tminus
                        elif V[i+1]>V[i-1]: # At local extremum, next V is higher
                            tangent = tplus*dVmax + tminus*dVmin
                        else: # At local extremum, previous V is higher
                            tangent = tplus*dVmin + tminus*dVmax
                        
                        # Normalize tangent
                        if np.linalg.norm(tangent) == 0: pass
                        else: tangent /= np.linalg.norm(tangent)
                    
                        # Find spring forces parallel to tangent
                        F_spring_parallel = NEB.k*( utils.dist(c,b) - utils.dist(b,a) ) * tangent
                        
                        energies[i] += 627.5 * 0.5*NEB.k*( utils.dist_squared(c,b) + utils.dist_squared(b,a) )
                    
                        # Find DFT forces perpendicular to tangent
                        real_force = np.array( [b.fx,b.fz,b.fz] )
                        F_real_perpendicular = real_force - np.dot(real_force, tangent)*tangent
                    
                        # Sum convergence check
                        NEB.convergence += b.fx**2 + b.fz**2 + b.fz**2
                        
                        # Set NEB forces
                        b.fx, b.fy, b.fz = F_spring_parallel + F_real_perpendicular
            # Get the RMF real force
            NEB.convergence = NEB.convergence**0.5

            # Set gradient
            NEB.gradient = []
            for state in NEB.states[1:-1]:
                for a in state:
                    NEB.gradient += [-a.fx, -a.fy, -a.fz] # Gradient of NEB.error

            # Calculate RMS Force
            RMS_force = (sum([a.fx**2+a.fy**2+a.fz**2 for state in states[1:-1] for a in state])/len([a for state in states[1:-1] for a in state]))**0.5
            NEB.RMS_force = RMS_force
                        
            # Print data
            V = V[:1] + [ (e-V[0])/0.001 for e in V[1:] ]
            if NEB.prv_RMS == None or NEB.prv_RMS > RMS_force:
                rms = utils.color_set(RMS_force,'GREEN')
            else:
                rms = utils.color_set(RMS_force,'RED')
            print NEB.step, '%7.5g +' % V[0], ('%5.1f '*len(V[1:])) % tuple(V[1:]), rms
            
            # Set error
            #NEB.error = max(V)
            NEB.error = RMS_force
            #NEB.error = 0.2*max(V) + 0.8*RMS_force
            #NEB.error = 0.5*sum([a.fx**2+a.fy**2+a.fz**2 for state in states[1:-1] for a in state]) # Max(V) # Sum(energies)

            NEB.prv_RMS = RMS_force
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

    def recenter(r):
        # Prevent rotation or translation
        coord_count = 0
        st = NEB.states
        for s in st[1:-1]:
            for a in s:
                a.x, a.y, a.z = r[coord_count], r[coord_count+1], r[coord_count+2]
                coord_count += 3
        utils.procrustes(st) #translate and rotate each frame to fit its neighbor
        coord_count = 0
        for s in st[1:-1]:
            for a in s:
                r[coord_count:coord_count+3] = [a.x, a.y, a.z]
                coord_count += 3

        return r

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

    def steepest_decent(f, r, fprime, alpha=0.1, maxiter=1000): #better, but tends to push error up eventually, especially towards endpoints.
        for step in range(maxiter+1):
            if NEB.convergence < NEB.convergence_criteria:
                print("\nConvergence achieved in %d iterations with %lg < %lg\n" % (NEB.step,NEB.convergence,NEB.convergence_criteria))
                sys.exit()
            gradient = -fprime(r)
            r += gradient*alpha
            if center:
                r = recenter(r)

    def quick_min_optimizer(f, r, nframes, fprime, dt=0.1, step_max=0.1, euler=False, viscosity=0.1, maxiter=1000): # dt = fs, step_max = angstroms, viscosity = 1/fs
        v = np.array([0.0 for x in r])
        acc = np.array([0.0 for x in r])

        masses = []
        for s in states[1:-1]:
            for a in s:
                m = units.elem_weight(a.element)
                masses += [m, m, m]
        masses = np.array(masses)

        for step in range(maxiter+1):
            if NEB.convergence < NEB.convergence_criteria:
                print("\nConvergence achieved in %d iterations with %lg < %lg\n" % (NEB.step,NEB.convergence,NEB.convergence_criteria))
                sys.exit()
            
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
            
            if center:
                r = recenter(r)

    def fire_optimizer(f, r, nframes, fprime, dt = 0.1, dtmax = 1.0, step_max = 0.2, maxiter=1000, 
                        Nmin = 5, finc = 1.1, fdec = 0.5, astart = 0.1, fa = 0.99, euler = True):

        v = np.array([0.0 for x in r])
        Nsteps = 0
        acc = astart

        for step in range(maxiter+1):
            if NEB.convergence < NEB.convergence_criteria:
                print("\nConvergence achieved in %d iterations with %lg < %lg\n" % (NEB.step,NEB.convergence,NEB.convergence_criteria))
                sys.exit()

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

            if center:
                r = recenter(r)

    def bfgs_optimizer(f, x0, fprime,
            alpha=0.1, beta=0.6, tau=1E-3, H_reset=True, gtol=1E-5,
            MIN_STEP=1E-10, MAX_STEP=0.2, reset=60, MAX_BACKTRACK=None,
            maxiter=1000, linesearch='backtrack', L2norm=True,
            BACKTRACK_EPS=1E-3, disp=0, callback=None):
        import numpy as np
        from math import (copysign, fsum)

        def vecnorm(x, ord=2):
            if ord == np.Inf:
                return np.amax(np.abs(x))
            elif ord == -np.Inf:
                return np.amin(np.abs(x))
            else:
                return np.sum(np.abs(x)**ord, axis=0)**(1.0 / ord)

        # Func Call Counter
        fcount = 0

        # Set maxiter if not set
        if maxiter is None:
            maxiter = 200 * len(x0)

        # Ensure coordinates are in the correct format
        x0 = np.asarray(x0).flatten()
        if x0.ndim == 0:
            x0.shape = (1,)

        # Initialize inv Hess and Identity matrix
        I = np.eye(len(x0), dtype=int)
        Hk = I

        # Get gradient and store your old func_max
        gfk = fprime(x0)
        if f is not None:
            old_fval = f(x0)
        fcount += 1
        # Hold coordinates for loop
        xk = x0

        # Criteria on gradient for continuing simulation
        if L2norm:
            def f_gnorm(x): 
                return np.sqrt(fsum([np.linalg.norm(g)**2 for g in x]))
            gnorm = f_gnorm(gfk)
        else:
            f_gnorm = vecnorm
            gnorm = f_gnorm(gfk, ord=norm)

        # Hold original values
        ALPHA_CONST = alpha
        BETA_CONST = beta
        RESET_CONST = reset

        # Get function to describe linesearch
        if linesearch is 'backtrack':
            if disp > 1:
                print("Using backtrack method with backtrack epsilon %lg." %
                    BACKTRACK_EPS)

            def check_backtrack(f1,f0,x,y):
                return (f1-f0)/(abs(f1)+abs(f0)) > BACKTRACK_EPS

        elif linesearch is 'armijo':
            if disp > 1:
                print("Using armijo method with tau %lg." % tau)

            def check_backtrack(f1,f0,pk,gk):
                return f1-f0 > tau*alpha*np.dot(gk,pk)

        else:
            print("ERROR - Linesearch '%s' is not acceptable. Use 'backtrack'\
                or 'armijo'")
            sys.exit()

        backtrack, loop_counter, warnflag = 0, 0, 0
        while (gnorm > gtol) and (loop_counter < maxiter):
            if MAX_BACKTRACK is not None and backtrack > MAX_BACKTRACK:
                warnflag = 2
                break
            if disp > 1:
                print("Step %d, " % loop_counter),

            # Get your step direction
            pk = -np.dot(Hk, gfk)

            # If we are doing unreasonably small step sizes, quit
            if abs(max(pk*alpha)) < MIN_STEP:
                if disp > 1:
                    print("Error - Step size unreasonable (%lg)" 
                                % abs(max(pk*alpha))),
                warnflag = 2
                break

            # If we have too large of a step size, set to max
            if max([abs(p*alpha) for p in pk]) > MAX_STEP:
                if disp > 1:
                    print("Warning - Setting step to max step size"
                                % np.linalg.norm(pk*alpha)),
                for i,p in enumerate(pk):
                    if abs(p*alpha) > MAX_STEP:
                        pk[i] = (MAX_STEP / alpha) * copysign(1, p)
                # As we are changing values manually, this is no longer
                # the BFGS(Hess) algorithm so reset the Inverse Hessian
                if H_reset:
                    Hk = I

            # Hold new parameters
            xkp1 = xk + alpha * pk

            # Recenter position
            if center:
                xkp1 = recenter(xkp1)

            # Recalculate sk to maintain the secant condition
            sk = xkp1 - xk

            # Get the new gradient
            gfkp1 = fprime(xkp1)
            # Check if max has increased
            if f is not None:
                fval = f(xkp1)
            fcount += 1

            if f is not None and check_backtrack(fval, old_fval, gfkp1, pk):
                # Step taken overstepped the minimum.  Lowering step size
                if disp > 1:
                    print("\tResetting System as %lg > %lg!"
                            % (fval, old_fval))
                    print("\talpha: %lg" % alpha),

                alpha *= np.float64(beta)

                if disp > 1:
                    print("-> %lg\n" % alpha)

                # Reset the Inverse Hessian if desired - This is recommended!
                if H_reset:
                    Hk = I
                backtrack += 1
                reset = RESET_CONST
                continue

            # This allows for the edge case in which after decreasing alpha, a situation arises
            # in which larger alphas are acceptable again. Thus, we reset to the original alpha
            elif reset is not None:
                reset -= 1
                if reset < 0 and alpha < ALPHA_CONST:
                    if disp > 1:
                        print("\tResetting Alpha, Beta, Reset and Inverse Hessian")
                    alpha, beta, reset = ALPHA_CONST, BETA_CONST, RESET_CONST
                    if H_reset:
                        Hk = I
                    continue
                elif reset < 0 and alpha >= ALPHA_CONST:
                    if disp > 1:
                        print("\tIncreasing step size: %lg ->" % alpha),
                    alpha /= beta
                    if disp > 1:
                        print("%lg,\t" % alpha),
            
            # Store new parameters, as it has passed the check
            # (fval < old_fval is True)
            xk = xkp1

            # Store new max value in old_max for future comparison
            if f is not None:
                old_fval = fval

            # Get difference in gradients for further calculations
            yk = gfkp1 - gfk
            # Store new gradient in old gradient
            gfk = gfkp1

            # Update the conditional check
            if L2norm:
                gnorm = f_gnorm(gfk)
            else:
                gnorm = f_gnorm(gfk, ord=norm)

            if disp > 1:
                print("gnorm %lg, fval %lg" % (gnorm, fval))

            # If callback is desired
            if callback is not None:
                callback(xk)

            # Increment the loop counter
            loop_counter += 1

            if L2norm:
                gnorm = f_gnorm(gfk)
            else:
                gnorm = f_gnorm(gfk, ord=norm)

            if (gnorm <= gtol):
                break

            try:  # this was handled in numeric, let it remaines for more safety
                rhok = 1.0 / (np.dot(yk, sk))
            except ZeroDivisionError:
                rhok = 1000.0
                if disp > 1:
                    print("Divide-by-zero encountered: rhok assumed large")
            if np.isinf(rhok):  # this is patch for np
                rhok = 1000.0
                if disp > 1:
                    print("Divide-by-zero encountered: rhok assumed large")
            # Run BFGS Update for the Inverse Hessian
            A1 = I - sk[:, np.newaxis] * yk[np.newaxis, :] * rhok
            A2 = I - yk[:, np.newaxis] * sk[np.newaxis, :] * rhok
            Hk = np.dot(A1, np.dot(Hk, A2)) + \
                 (rhok * sk[:, np.newaxis] * sk[np.newaxis, :])

        if f is not None:
            fval = old_fval
        else:
            fval = float('NaN')

        if np.isnan(fval):
            # This can happen if the first call to f returned NaN;
            # the loop is then never entered.
            warnflag = 2

        if warnflag == 2:
            if disp == 1:
                print("Warning: Loss of precision.")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % fcount)

        elif loop_counter >= maxiter:
            warnflag = 1
            if disp == 1:
                print("Warning: Maximum Iteration was exceeded.")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % fcount)
        else:
            if disp == 1:
                print("Success!")
                print("         Current function value: %f" % fval)
                print("         Iterations: %d" % loop_counter)
                print("         Function evaluations: %d" % fcount)

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

    n = NEB(name, states, theory, extra_section, k, fit_rigid)

    # Output for user
    if opt == 'BROYDEN_ROOT':
        scipy.optimize.broyden1(NEB.get_gradient, np.array(NEB.coords_start), alpha=float(alpha), verbose=(disp != 0))
    elif opt == 'QM':
        quick_min_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.nframes, 
            fprime=NEB.get_gradient, dt=dt, viscosity=viscosity, step_max=step_max, euler=euler, maxiter=maxiter)
    elif opt == 'FIRE':
        fire_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.nframes, 
            fprime=NEB.get_gradient, dt=dt, dtmax=dtmax, step_max=step_max,
            Nmin=Nmin, finc=finc, fdec=fdec, astart=astart, fa=fa, euler=euler, maxiter=maxiter)
    elif opt == 'BFGS':
        bfgs_optimizer(NEB.get_error, np.array(NEB.coords_start), fprime=NEB.get_gradient,
            gtol=float(gtol), maxiter=int(maxiter),
            alpha=float(alpha), beta=float(beta), tau=float(tau), H_reset=H_reset,
            MIN_STEP=float(step_min), MAX_STEP=float(step_max), reset=reset, MAX_BACKTRACK=bt_max, 
            linesearch=linesearch, L2norm=L2norm, BACKTRACK_EPS=bt_eps, disp=disp
            )
    elif opt == 'BFGS2':
        from scipy.optimize.bfgsh import fmin_bfgs_h
        fmin_bfgs_h(NEB.get_error, np.array(NEB.coords_start), fprime=NEB.get_gradient,
            alpha=float(alpha), beta=float(beta), tau=float(tau), H_reset=H_reset,
            MIN_STEP=float(step_min), MAX_STEP=float(step_max), reset=reset, MAX_BACKTRACK=bt_max, 
            linesearch=linesearch, L2norm=L2norm, BACKTRACK_EPS=bt_eps, disp=(disp>0)
            )
    elif opt == 'SD':
        steepest_decent(NEB.get_error, np.array(NEB.coords_start), fprime=NEB.get_gradient, alpha=alpha, maxiter=maxiter)
    else:
        print("\nERROR - %s optimizations method does not exist! Choose from the following:" % str(opt))
        print("\t1. BFGS")
        print("\t2. BFGS2")
        print("\t3. QM")
        print("\t4. SD")
        print("\t5. FIRE")
        print("\t6. BROYDEN_ROOT\n")
        sys.exit()
