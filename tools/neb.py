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
        fit_rigid=True, center=True, k=0.1837,
        DFT='g09', opt='BFGS', scipy_test=False, gtol=1e-5,
        alpha=0.3, beta=0.4, H_reset=True, dt = 0.5, euler=True,
        force=True, mem=25, blurb=None, initial_guess=None): 
    
    # If using test code, import path so we import correct scipy.optimize.
    if scipy_test or opt=='BFGS2': sys.path.insert(1,'/fs/home/hch54/scipy_mod/scipy/')
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
    elif opt in ['BFGS','BFGS2']:
        print("\nRunning neb with optimization method %s" % str(opt))
        print("\talpha = %lg, beta = %lg, H_reset = %s" % (alpha, beta, str(H_reset)))
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
                if NEB.step>0:
                    if (i==0 or i==len(NEB.states)-1): # Only use DFT on endpoints on first step, because they don't change
                        continue
                    if DFT=='g09':
                        guess = ' Guess=Read'
                    elif DFT=='orca':
                        guess = ' MOREAD'
                        tmp = '%%moinp "../%s-%d-%d/%s-%d-%d.orca.gbw"' % (NEB.name,NEB.step-1,i,NEB.name,NEB.step-1,i)
                        extra_section = tmp + NEB.extra_section.strip()
                else:
                    if initial_guess:
                        if DFT=='g09':
                            guess = ' Guess=Read'
                        elif DFT=='orca':
                            guess = ' MOREAD'
                            tmp = '%%moinp "../%s/%s.orca.gbw\n"' % (initial_guess,initial_guess)
                            extra_section = tmp + NEB.extra_section.strip()
                    else:
                        guess = '' #no previous guess for first step
                        extra_section =  NEB.extra_section.strip()
                if DFT=='g09':
                    if extra_section != '':
                        extra_section = extra_section.strip() + '\n\n'
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
            
            # Set error
            #NEB.error = units.convert_energy('Ha','kcal/mol',max(energies)-energies[0])
            #NEB.error = 0.5*sum([a.fx**2+a.fy**2+a.fz**2 for state in states[1:-1] for a in state]) # Max(V) # Sum(energies)
            NEB.error = RMS_force
            
            # Print data
            V = V[:1] + [ (e-V[0])/0.001 for e in V[1:] ]
            if NEB.prv_RMS == None or NEB.prv_RMS > RMS_force:
                rms = utils.color_set(RMS_force,'GREEN')
            else:
                rms = utils.color_set(RMS_force,'RED')
            print NEB.step, '%7.5g +' % V[0], ('%5.1f '*len(V[1:])) % tuple(V[1:]), rms
            
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
            v1, v2: numpy vectors
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

    def steepest_decent(f, r, fprime, alpha=0.1): #better, but tends to push error up eventually, especially towards endpoints.
        print("Running with alpha = %lg" % alpha)
        for step in range(1000):
            if NEB.convergence < NEB.convergence_criteria:
                print("\nConvergence achieved in %d iterations with %lg < %lg\n" % (NEB.step,NEB.convergence,NEB.convergence_criteria))
                sys.exit()
            print("%d. Real RMS: %lg," % (NEB.step,NEB.convergence)),
            gradient = -fprime(r)
            r += gradient*alpha
            if center:
                r = recenter(r)

    def quick_min_optimizer(f, r, nframes, fprime, dt=0.5, max_dist=0.1, euler=False, viscosity=0.1): # dt = fs, max_dist = angstroms, viscosity = 1/fs
        v = np.array([0.0 for x in r])
        acc = np.array([0.0 for x in r])

        masses = []
        for s in states[1:-1]:
            for a in s:
                m = units.elem_weight(a.element)
                masses += [m, m, m]
        masses = np.array(masses)

        for step in range(1000):
            if NEB.convergence < NEB.convergence_criteria:
                print("\nConvergence achieved in %d iterations with %lg < %lg\n" % (NEB.step,NEB.convergence,NEB.convergence_criteria))
                sys.exit()
            print("%d. Real RMS: %lg," % (NEB.step,NEB.convergence)),

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
                if speed*dt > max_dist:
                    max_speed = max_dist/dt
                    v[low:high] *= max_speed / speed
                
            if euler:
                #make Euler step
                v += dt * forces

                #limit distance moved
                #for i in range(len(v)):
                #   if v[i]*dt > max_dist: v[i] = max_dist/dt
                #   if v[i]*dt <-max_dist: v[i] =-max_dist/dt

                for i in range(1,nframes-1):
                    low = (i-1)*natoms*3
                    high = i*natoms*3
                    speed = np.linalg.norm(v[low:high])
                    if speed*dt > max_dist:
                        max_speed = max_dist/dt
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
                #   if dx[i] > max_dist: dx[i] = max_dist
                #   if dx[i] <-max_dist: dx[i] =-max_dist
                
                r_new = r + dx
                v_new = v + (acc + a_new)*0.5 * dt
                r = r_new
                v = v_new
                acc = a_new
            
            if center:
                r = recenter(r)

    def fire_optimizer(f, r, nframes, fprime, dt = 0.05, dtmax = 1.0, max_dist = 0.2, 
                        Nmin = 5, finc = 1.1, fdec = 0.5, astart = 0.1, fa = 0.99, euler = True):

        v = np.array([0.0 for x in r])
        Nsteps = 0
        acc = astart

        for step in range(1000):
            if NEB.convergence < NEB.convergence_criteria:
                print("\nConvergence achieved in %d iterations with %lg < %lg\n" % (NEB.step,NEB.convergence,NEB.convergence_criteria))
                sys.exit()
            print("%d. Real RMS: %lg," % (NEB.step,NEB.convergence)),
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
                    if v[i]*dt > max_dist: v[i] = max_dist/dt
                    if v[i]*dt <-max_dist: v[i] =-max_dist/dt

                #move atoms
                r += v * dt

            if center:
                r = recenter(r)

    def bfgs_optimizer(f, r, fprime, alpha=0.3, beta=0.7, H_reset=True, gtol=1e-5):
        # BFGS optimizer adapted from scipy.optimize._minmize_bfgs
        import numpy as np
        def vecnorm(x, ord=2):
            if ord == np.Inf:
                return np.amax(np.abs(x))
            elif ord == -np.Inf:
                return np.amin(np.abs(x))
            else:
                return np.sum(np.abs(x)**ord, axis=0)**(1.0 / ord)

        if beta > 1:
            print("Warning - Unreasonable Beta (must be less than or equal to 1). Setting to 1.\n")
            beta = 1.0

        # Get x0 as a flat array
        x0 = np.asarray(r).flatten()
        if x0.ndim == 0:
            x0.shape = (1,)
        
        maxiter = len(x0) * 200

        g0 = fprime(x0)

        loop_counter,N = 0,len(x0)
        I = np.eye(N, dtype=int)

        # Initialize inv Hess as Identity matrix
        Hk = I

        # Store your old energy
        E_old = f(x0)

        xk = x0

        # Criteria on gradient for continuing simulation
        norm = np.Inf
        gnorm = vecnorm(g0, ord=norm)
        
        # Main Loop:
        backtrack_counter = 0
        while (gnorm > gtol) and (loop_counter < maxiter):
            if backtrack_counter >= 10: break
            print("%d. Real RMS: %lg," % (NEB.step,NEB.convergence)),
            # Get your step direction
            pk = -np.dot(Hk, g0)

            # If we are doing unreasonably small step sizes, quit
            if np.linalg.norm(pk*alpha) < 1E-7:
                print("Error - Step size unreasonable (%lg Angstroms)" % np.linalg.norm(pk*alpha))
                sys.exit()

            # Hold new position
            xkp1 = xk + alpha * pk

            # Recalculate sk to maintain the secant condition
            sk = xkp1 - xk

            # Recenter position
            if center:
                xkp1 = recenter(xkp1)

            # Get the new gradient
            gfkp1 = fprime(xkp1)

            # Check if max energy has increased
            E_new = f(xkp1)
            if E_new > E_old:
                # Step taken overstepped the minimum.  Lowering step size
                print("Resetting System as %lg > %lg!" % (E_new, E_old))
                print("\talpha: %lg" % alpha),
                alpha *= float(beta)
                print("-> %lg" % alpha)
                print("\tmag(sk) = %lg" % vecnorm(sk, ord=norm))
                if 'yk' in locals():
                    print("\t<yk|sk> = %lg\n" % (np.dot(yk, sk)))
                # NOTE! Maybe add a scalar for lowering the mag of Hk
                if H_reset: Hk = I
                backtrack_counter += 1
                continue
            
            # Store new position, as it has passed the check (E_new < E_old is True)
            xk = xkp1
            
            # Store new energy in old energy for future comparison
            E_old = E_new

            # Get difference in gradients for further calculations
            yk = gfkp1 - g0
            # Store new gradient in old gradient
            g0 = gfkp1

            try:  # this was handled in numeric, let it remaines for more safety
                rhok = 1.0 / (np.dot(yk, sk))
            except ZeroDivisionError:
                rhok = 1000.0
                print("Divide-by-zero encountered: rhok assumed large")
            if np.isinf(rhok):  # this is patch for numpy
                rhok = 1000.0
                print("Divide-by-zero encountered: rhok assumed large")

            # Run BFGS Update for the Inverse Hessian
            A1 = I - sk[:, np.newaxis] * yk[np.newaxis, :] * rhok
            A2 = I - yk[:, np.newaxis] * sk[np.newaxis, :] * rhok
            Hk = np.dot(A1, np.dot(Hk, A2)) + \
                 (rhok * sk[:, np.newaxis] * sk[np.newaxis, :])

            # Update the conditional check
            gnorm = vecnorm(g0, ord=norm)

            # Increment the loop counter
            loop_counter += 1

        if NEB.convergence < NEB.convergence_criteria:
            print("\nConvergence achieved in %d iterations with %lg < %lg\n" % (NEB.step,NEB.convergence,NEB.convergence_criteria))
            sys.exit()

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
        scipy.optimize.broyden1(NEB.get_gradient, np.array(NEB.coords_start), alpha=float(alpha), verbose=True)
    elif opt == 'QM':
        quick_min_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.nframes, 
            fprime=NEB.get_gradient, dt=dt, max_dist=0.01, euler=euler)
    elif opt == 'FIRE':
        fire_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.nframes, 
            fprime=NEB.get_gradient, dt = dt, dtmax = 1.0, max_dist = 0.2,
            Nmin = 5, finc = 1.1, fdec = 0.5, astart = 0.1, fa = 0.99, euler=euler)
    elif opt == 'BFGS':
        bfgs_optimizer(NEB.get_error, np.array(NEB.coords_start), fprime=NEB.get_gradient, alpha=float(alpha), beta=float(beta), gtol=gtol, H_reset=H_reset)
    elif opt == 'BFGS2':
        from scipy.optimize.bfgsh import fmin_bfgs_h
        fmin_bfgs_h(NEB.get_error, np.array(NEB.coords_start), fprime=NEB.get_gradient, alpha=float(alpha), beta=float(beta), gtol=gtol, H_reset=H_reset)
    elif opt == 'SD':
        steepest_decent(NEB.get_error, np.array(NEB.coords_start), fprime=NEB.get_gradient, alpha=alpha)
    else:
        print("\nERROR - %s optimizations method does not exist! Choose from the following:" % str(opt))
        print("\t1. BFGS")
        print("\t2. BFGS2")
        print("\t3. QM")
        print("\t4. SD")
        print("\t5. FIRE")
        print("\t6. BROYDEN_ROOT\n")
        sys.exit()
