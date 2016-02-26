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

import numpy as np
import math

def euler2mat(z=0, y=0, x=0):
    ''' Return matrix for rotations around z, y and x axes

    Uses the z, then y, then x convention above

    Parameters
    ----------
    z : scalar
       Rotation angle in radians around z-axis (performed first)
    y : scalar
       Rotation angle in radians around y-axis
    x : scalar
       Rotation angle in radians around x-axis (performed last)

    Returns
    -------
    M : array shape (3,3)
       Rotation matrix giving same rotation as for given angles

    Examples
    --------
    >>> zrot = 1.3 # radians
    >>> yrot = -0.1
    >>> xrot = 0.2
    >>> M = euler2mat(zrot, yrot, xrot)
    >>> M.shape == (3, 3)
    True

    The output rotation matrix is equal to the composition of the
    individual rotations

    >>> M1 = euler2mat(zrot)
    >>> M2 = euler2mat(0, yrot)
    >>> M3 = euler2mat(0, 0, xrot)
    >>> composed_M = np.dot(M3, np.dot(M2, M1))
    >>> np.allclose(M, composed_M)
    True

    You can specify rotations by named arguments

    >>> np.all(M3 == euler2mat(x=xrot))
    True

    When applying M to a vector, the vector should column vector to the
    right of M.  If the right hand side is a 2D array rather than a
    vector, then each column of the 2D array represents a vector.

    >>> vec = np.array([1, 0, 0]).reshape((3,1))
    >>> v2 = np.dot(M, vec)
    >>> vecs = np.array([[1, 0, 0],[0, 1, 0]]).T # giving 3x2 array
    >>> vecs2 = np.dot(M, vecs)

    Rotations are counter-clockwise.

    >>> zred = np.dot(euler2mat(z=np.pi/2), np.eye(3))
    >>> np.allclose(zred, [[0, -1, 0],[1, 0, 0], [0, 0, 1]])
    True
    >>> yred = np.dot(euler2mat(y=np.pi/2), np.eye(3))
    >>> np.allclose(yred, [[0, 0, 1],[0, 1, 0], [-1, 0, 0]])
    True
    >>> xred = np.dot(euler2mat(x=np.pi/2), np.eye(3))
    >>> np.allclose(xred, [[1, 0, 0],[0, 0, -1], [0, 1, 0]])
    True

    Notes
    -----
    The direction of rotation is given by the right-hand rule (orient
    the thumb of the right hand along the axis around which the rotation
    occurs, with the end of the thumb at the positive end of the axis;
    curl your fingers; the direction your fingers curl is the direction
    of rotation).  Therefore, the rotations are counterclockwise if
    looking along the axis of rotation from positive to negative.
    '''
    Ms = []
    if z:
        cosz = math.cos(z)
        sinz = math.sin(z)
        Ms.append(np.array(
                [[cosz, -sinz, 0],
                 [sinz, cosz, 0],
                 [0, 0, 1]]))
    if y:
        cosy = math.cos(y)
        siny = math.sin(y)
        Ms.append(np.array(
                [[cosy, 0, siny],
                 [0, 1, 0],
                 [-siny, 0, cosy]]))
    if x:
        cosx = math.cos(x)
        sinx = math.sin(x)
        Ms.append(np.array(
                [[1, 0, 0],
                 [0, cosx, -sinx],
                 [0, sinx, cosx]]))
    if Ms:
        return reduce(np.dot, Ms[::-1])
    return np.eye(3)



def sum_NEB_springs(translations_and_rotations, V, NEB, energies, fake_states=None, pointer=None):
    N_modified = len(NEB.states)-1 #skip first image
    #unpack translations
    translations = []
    for i in range(0, N_modified*3, 3): #three degrees of freedom per frame
        translations.append( translations_and_rotations[i:i+3] )
    #unpack rotations
    full_rotation = []
    for i in range(N_modified*3, N_modified*3*2, 3):
        euler_angles = translations_and_rotations[i:i+3]
        full_rotation.append( euler2mat(*euler_angles) )
        
    #make rotation matrices from 3 angles
    
    if not fake_states:
        fake_states = copy.deepcopy(NEB.states)
    
    #start_rotation = utils.procrustes(fake_states, append_in_loop=False)
    #
    
    #for i in range(len(full_rotation)):
    #    full_rotation[i] = np.dot( full_rotation[i], start_rotation[i] )
    
    for i,state in enumerate(fake_states[1:]):
        for a in state:
            a.x += translations[i][0]
            a.y += translations[i][1]
            a.z += translations[i][2]
            a.x,a.y,a.z = np.dot( full_rotation[i], (a.x,a.y,a.z) )
    
    sum_springs = 0.0
    for i,state in enumerate(NEB.states):
        if i==0 or i==len(NEB.states)-1: continue # Don't change first or last state
        for j,bb in enumerate(state):
            if j in range(len(NEB.states[0])):
                F_spring_parallel, F_real_perpendicular = get_NEB_forces(fake_states, i, j, V, full_rotation, NEB, energies, bb, False)
                sum_springs += np.linalg.norm(F_spring_parallel)
    
    if pointer is not None:
        pointer.append(full_rotation)
    
    return sum_springs# + 1e3*sum([x**2 for x in translations_and_rotations])





def get_NEB_forces(fake_states, i, j, V, full_rotation, NEB, energies, bb, frigid):
    import numpy as np
    a,b,c = [fake_states[i+d][j] for d in [-1,0,1]]

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
    if np.linalg.norm(tangent) != 0:
        tangent /= np.linalg.norm(tangent)

    # Find spring forces parallel to tangent
    F_spring_parallel = NEB.k*( utils.dist(c,b) - utils.dist(b,a) ) * tangent

    E_hartree = 0.5*NEB.k*( utils.dist_squared(c,b) + utils.dist_squared(b,a) )
    energies[i] += units.convert_energy('Ha','kcal/mol', E_hartree)

    # Find DFT forces perpendicular to tangent
    real_force = np.array( [bb.fx,bb.fy,bb.fz] )
    F_real_perpendicular = real_force - np.dot(real_force, tangent)*tangent

    # Sum convergence check
    NEB.convergence += b.fx**2 + b.fy**2 + b.fz**2
    
    if frigid:
        R = full_rotation[i-1]
        R_inv = np.linalg.inv(R)
        F_spring_parallel = np.dot( R_inv, F_spring_parallel )
    
    return F_spring_parallel, F_real_perpendicular





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
        guess = ' MOREAD'
        tmp = '%%moinp "../%s-%d-%d/%s-%d-%d.orca.gbw"' % (NEB.name,NEB.step-1,i,NEB.name,NEB.step-1,i)
        extra_section = tmp + extra_section.strip()
    else:
        if initial_guess:
            guess = ' MOREAD'
            tmp = '%%moinp "../%s/%s.orca.gbw\n"' % (initial_guess,initial_guess)
            extra_section = tmp + extra_section.strip()
        else:
            guess = '' #no previous guess for first step
    return orca.job('%s-%d-%d'%(NEB.name,NEB.step,i), NEB.theory+guess, state, extra_section=extra_section, grad=True, procs=procs, queue=queue)


def orca_results(NEB, step_to_use, i, state):
    new_atoms, new_energy = orca.engrad_read('%s-%d-%d' % (NEB.name, step_to_use, i),force='Ha/Ang',pos='Ang')
    
    for a,b in zip(state, new_atoms):
        a.fx,a.fy,a.fz = b.fx,b.fy,b.fz
    
    return new_energy, new_atoms


def neb(name, states, theory, extra_section='', spring_atoms=None, procs=1, queue=None,
        disp=0, k=0.1837, frigid=True, start_job=None, get_results=None,
        DFT='orca', opt='LBFGS', gtol=1e-3, maxiter=1000,
        alpha=0.1, beta=0.5, tau=1E-3, reset=10, H_reset=True, Nmax=20,
        viscosity=0.1, dtmax=1.0, Nmin=5, finc=1.1, fdec=0.5, astart=0.1, fa=0.99,
        step_min=1E-8, step_max=0.2, bt_max=None, linesearch='backtrack', L2norm=True, bt_eps=1E-3,
        dt = 0.3, euler=True, force=True, mem=25, blurb=None, initial_guess=None): 
    
    import scipy.optimize
    import numpy as np

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
        def __init__(self, name, states, theory, extra_section='', k=0.1837):
            NEB.name = name
            NEB.states = states
            NEB.theory = theory
            NEB.extra_section = extra_section
            NEB.k = k
            NEB.prv_RMS = None
            NEB.convergence_criteria = gtol
            NEB.convergence = float('inf')
            NEB.nframes = len(states)
            NEB.RMS_force = float('inf')

            utils.procrustes(NEB.states) # In all cases, optimize the path via rigid rotations first
            
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
                    pass #No need to calculate anything for first and last states after the first step
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
            # Reset convergence check
            NEB.convergence = 0

            fake_states = copy.deepcopy(NEB.states)
            fake_states2 = copy.deepcopy(NEB.states)
            #full_rotation = utils.procrustes(fake_states2, append_in_loop=False)
            
            #print sum_NEB_springs([0.0 for x in range( (len(NEB.states)-1)*6 )], V, NEB, energies)
            solution = scipy.optimize.minimize(sum_NEB_springs, [0.0 for x in range( (len(NEB.states)-1)*6 )], args=(V, NEB, energies), options={'disp':False}, method='Nelder-Mead')
            #print solution.x
            pointer = []
            sum_NEB_springs(solution.x, V, NEB, energies, fake_states2, pointer)
            full_rotation = pointer[0]
            
            #if frigid:
            #    full_rotation = utils.procrustes(fake_states, append_in_loop=False)
            #else:
            #    full_rotation = None
            
            sum_spring, sum_spring2 = 0.0, 0.0
            
            # Add spring forces to atoms
            for i,state in enumerate(NEB.states):
                if i==0 or i==len(NEB.states)-1: continue # Don't change first or last state
                for j,bb in enumerate(state):
                    if j in spring_atoms:
                        F_spring_parallel, F_real_perpendicular = get_NEB_forces(fake_states, i, j, V, full_rotation, NEB, energies, bb, False)
                        F_spring_parallel2, F_real_perpendicular2 = get_NEB_forces(fake_states2, i, j, V, full_rotation, NEB, copy.deepcopy(energies), bb, True)
                        
                        sum_spring += np.linalg.norm(F_spring_parallel)
                        sum_spring2 += np.linalg.norm(F_spring_parallel2)
                        
                        # Set NEB forces
                        if frigid:
                            bb.fx, bb.fy, bb.fz = F_spring_parallel2 + F_real_perpendicular2
                        else:
                            bb.fx, bb.fy, bb.fz = F_spring_parallel + F_real_perpendicular
                        
            print sum_spring, sum_spring2
            #remove net translation forces from the gradient
            
            if frigid != None:
                for state in NEB.states[1:-1]:
                    net_force = np.zeros(3)
                    for a in state:
                        net_force += (a.fx, a.fy, a.fz)
                    for a in state:
                        a.fx -= net_force[0] / len(state)
                        a.fy -= net_force[1] / len(state)
                        a.fz -= net_force[2] / len(state)
            
            # Get the RMF real force
            NEB.convergence = NEB.convergence**0.5

            # Set gradient
            NEB.gradient = []
            for state in NEB.states[1:-1]:
                for a in state:
                    NEB.gradient += [-a.fx, -a.fy, -a.fz] # Gradient of NEB.error

            # Calculate RMS Force
            force_mags = [(a.fx**2+a.fy**2+a.fz**2)**0.5 for state in states[1:-1] for a in state]
            RMS_force = (sum([f**2 for f in force_mags])/float(len(force_mags)))**0.5
            NEB.RMS_force = RMS_force
                        
            # Print data
            V = V[:1] + [ (e-V[0])/0.001 for e in V[1:] ]
            if NEB.prv_RMS == None or NEB.prv_RMS > RMS_force:
                rms = utils.color_set(RMS_force,'GREEN')
            else:
                rms = utils.color_set(RMS_force,'RED')
            print NEB.step, '%7.5g +' % V[0], ('%5.1f '*len(V[1:])) % tuple(V[1:]), rms, 'Ha/Ang'
            
            if NEB.prv_RMS is None: NEB.prv_RMS = RMS_force
            NEB.prv_RMS = min(RMS_force, NEB.prv_RMS)

            # Set error
            #NEB.error = max(V)
            NEB.error = RMS_force
            #NEB.error = max(force_mags)
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
        from scipy.linalg import block_diag
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

    def steepest_descent(f, r, fprime, alpha=0.05, maxiter=1000, gtol=1E-3):
        step = 0
        while (NEB.RMS_force > gtol) and (step < maxiter):
            if NEB.convergence < NEB.convergence_criteria:
                print("\nConvergence achieved in %d iterations with %lg < %lg\n" % (NEB.step,NEB.convergence,NEB.convergence_criteria))
                sys.exit()
            forces = -fprime(r)
            r += forces*alpha
            #r = align_coordinates(r)
            step += 1

    def quick_min_optimizer(f, r, nframes, fprime, dt=0.1, step_max=0.1, euler=False, viscosity=0.1, maxiter=1000, gtol=1E-3): # dt = fs, step_max = angstroms, viscosity = 1/fs
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
            
            r = align_coordinates(r)

            step += 1

    def fire_optimizer(f, r, nframes, fprime, dt = 0.1, dtmax = 1.0, step_max = 0.2, maxiter=1000, gtol=1E-3,
                        Nmin = 5, finc = 1.1, fdec = 0.5, astart = 0.1, fa = 0.99, euler = True):

        v = np.array([0.0 for x in r])
        Nsteps = 0
        acc = astart

        step = 0
        while (NEB.RMS_force > gtol) and (step < maxiter):
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

            r = align_coordinates(r)

            step += 1


    def lbfgs_optimizer(target_function, initial_coordinates, target_gradient,
            step_size=0.1, step_size_adjustment=0.5, armijio_line_search_factor=1E-4, reset_when_in_trouble=True, linesearch='armijo',
            gradient_tolerance=1E-3, max_iterations=1000, reset_step_size=reset,
            MAX_STEP=0.2, fit_rigid=True,
            display=0, callback=None, max_steps_remembered=Nmax):
        import numpy as np

        # These are values deemed good for DFT NEB and removed from parameter space for simplification
        MAX_BACKTRACK=None
        MIN_STEP=1E-8
        BACKTRACK_EPS=1E-3

        if display > 2:
            print("\nValues in bfgs_optimize code:")
            print("\tstep size = %lg, step size adjustment = %lg, reset_when_in_trouble = %s" % (step_size, step_size_adjustment, str(reset_when_in_trouble)))
            print("\tgtol = %lg, max_iterations = %d, MAX_STEP = %lg" % (gradient_tolerance, max_iterations, MAX_STEP))
            print("\t-------------------")
            print("\tMAX_BACKTRACK = %s, reset_step_size = %s, MIN_STEP = %lg, BACKTRACK_EPS = %lg" % (str(MAX_BACKTRACK), str(reset_step_size), MIN_STEP, BACKTRACK_EPS))
            print("\tfrigid = %s\n" % str(fit_rigid))
        def vecnorm(x, ord=2):
            if ord == np.Inf:
                return np.amax(np.abs(x))
            elif ord == -np.Inf:
                return np.amin(np.abs(x))
            else:
                return np.sum(np.abs(x)**ord, axis=0)**(1.0 / ord)

        function_call_counter = 0

        # Set max_iterations if not set
        if max_iterations is None:
            max_iterations = 200 * len(initial_coordinates)

        # Ensure coordinates are in the correct format
        initial_coordinates = np.asarray(initial_coordinates).flatten()
        if initial_coordinates.ndim == 0:
            initial_coordinates.shape = (1,)
        current_coordinates = align_coordinates(initial_coordinates)

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
            if MAX_BACKTRACK is not None and backtrack > MAX_BACKTRACK:
                warnflag = 2
                break

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
        import numpy as np

        # These are values deemed good for DFT NEB and removed from parameter space for simplification
        MAX_BACKTRACK=None
        MIN_STEP=1E-8
        BACKTRACK_EPS=1E-3

        if display > 2:
            print("\nValues in bfgs_optimize code:")
            print("\tstep size = %lg, step size adjustment = %lg, reset_when_in_trouble = %s" % (step_size, step_size_adjustment, str(reset_when_in_trouble)))
            print("\tgtol = %lg, max_iterations = %d, MAX_STEP = %lg" % (gradient_tolerance, max_iterations, MAX_STEP))
            print("\t-------------------")
            print("\tMAX_BACKTRACK = %s, reset_step_size = %s, MIN_STEP = %lg, BACKTRACK_EPS = %lg" % (str(MAX_BACKTRACK), str(reset_step_size), MIN_STEP, BACKTRACK_EPS))
            print("\tfrigid = %s\n" % str(fit_rigid))
        def vecnorm(x, ord=2):
            if ord == np.Inf:
                return np.amax(np.abs(x))
            elif ord == -np.Inf:
                return np.amin(np.abs(x))
            else:
                return np.sum(np.abs(x)**ord, axis=0)**(1.0 / ord)

        function_call_counter = 0

        # Set max_iterations if not set
        if max_iterations is None:
            max_iterations = 200 * len(initial_coordinates)

        # Ensure coordinates are in the correct format
        initial_coordinates = np.asarray(initial_coordinates).flatten()
        if initial_coordinates.ndim == 0:
            initial_coordinates.shape = (1,)
        current_coordinates = align_coordinates(initial_coordinates)

        # Initialize inv Hess and Identity matrix
        I = np.eye(len(current_coordinates), dtype=int)
        current_Hessian = I

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
            if MAX_BACKTRACK is not None and backtrack > MAX_BACKTRACK:
                warnflag = 2
                break

            if display > 1:
                print("Trace of current_Hessian = %lg" % float(np.matrix.trace(current_Hessian)))
                print("Step %d, " % loop_counter),

            # Get your step direction
            step_direction = -np.dot(current_Hessian, current_gradient)
            # Renorm to remove the effect of current_Hessian not being unit
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
                    current_Hessian = I

            # Hold new parameters
            new_coordinates = current_coordinates + step_size * step_direction

            if fit_rigid:
                new_coordinates, C, current_Hessian_tmp = align_coordinates(new_coordinates, [current_gradient, current_coordinates], current_Hessian)
                current_gradient_tmp, current_coordinates_tmp = C

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
                # inverse hessian corects itself, it might not be important to do this.
                if reset_when_in_trouble:
                    current_Hessian = I
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
                        current_Hessian = I
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
            if target_function is not None:
                old_fval = fval

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

    n = NEB(name, states, theory, extra_section, k)

    # Output for user
    if opt == 'BROYDEN_ROOT':
        scipy.optimize.broyden1(NEB.get_gradient, np.array(NEB.coords_start), alpha=float(alpha), verbose=(disp != 0))
    elif opt == 'QM':
        quick_min_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.nframes, 
            fprime=NEB.get_gradient, dt=dt, viscosity=viscosity, step_max=step_max, euler=euler, maxiter=maxiter, gtol=gtol)
    elif opt == 'FIRE':
        fire_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.nframes, 
            fprime=NEB.get_gradient, dt=dt, dtmax=dtmax, step_max=step_max,
            Nmin=Nmin, finc=finc, fdec=fdec, astart=astart, fa=fa, euler=euler, maxiter=maxiter, gtol=gtol)
    elif opt == 'BFGS':
        bfgs_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.get_gradient,
            step_size=float(alpha), step_size_adjustment=float(beta), armijio_line_search_factor=float(tau), reset_when_in_trouble=H_reset,
            gradient_tolerance=float(gtol), max_iterations=int(maxiter), fit_rigid=frigid, linesearch=linesearch, reset_step_size=reset,
            MAX_STEP=float(step_max), display=disp
            )
    elif opt == 'LBFGS':
        lbfgs_optimizer(NEB.get_error, np.array(NEB.coords_start), NEB.get_gradient,
            step_size=float(alpha), step_size_adjustment=float(beta), armijio_line_search_factor=float(tau), reset_when_in_trouble=H_reset,
            gradient_tolerance=float(gtol), max_iterations=int(maxiter), fit_rigid=frigid, linesearch=linesearch, reset_step_size=reset,
            MAX_STEP=float(step_max), max_steps_remembered=Nmax, display=disp
            )
    elif opt == 'SD':
        steepest_descent(NEB.get_error, np.array(NEB.coords_start), fprime=NEB.get_gradient, alpha=alpha, maxiter=maxiter, gtol=gtol)
    else:
        print("\nERROR - %s optimizations method does not exist! Choose from the following:" % str(opt))
        print("\t1. BFGS")
        print("\t2. LBFGS")
        print("\t3. QM")
        print("\t4. SD")
        print("\t5. FIRE")
        print("\t6. BROYDEN_ROOT\n")
        sys.exit()
