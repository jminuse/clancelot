##################
# A BFGS optimizer
##################
#
# Required Values:
#   initial_parameters - Flat Array of floats - A set of initial parameters you wish to optimize.
#   gradient_function - Function - A function that take a list of floats (parameters) and returns a gradient to guide optimization.
#
# Optional Values:
#   target_function - Function - A function that, given a list of floats (parameters), returns a single value whose purpose is to verify
#                                that the parameters are minimizing appropriately.
#   dimension - int - 
#
#   step_size - float - A value specifying the initial step size to be used for .........
#   step_size_adjustment - float - 
#   max_step_size - float - 
#
#   fit_rigid - bool - 
#
#   armijo_line_search_factor - float - 
#   linesearch - str - 
#
#   max_gradient - float - 
#   rms_gradient - float - 
#   max_iterations - int - 
#
#   reset_step_size - int - 
#   reset_when_in_trouble - bool - 
#   
#   display - int - 
#   
#   callback - function - 
#
##################

from numpy import asarray, eye, sqrt, dot, float64, newaxis, isinf
from copy import deepcopy

def bfgs_optimizer(initial_parameters, gradient_function, target_function=None, dimension=3, start_hess=None,
                   step_size=0.1, step_size_adjustment=0.5, armijio_line_search_factor=1E-4, reset_when_in_trouble=True, linesearch='armijo',
                   max_gradient=1E-2, rms_gradient=1E-3, max_iterations=1000, reset_step_size=10,
                   max_step_size=0.2, fit_rigid=True, display=0, callback=None):

    # Wrap functions together
    def target(parameters, function_call_counter):
        function_call_counter += 1
        if target_function is not None:
            return target_function(parameters), gradient_function(parameters), function_call_counter
        else:
            return None, gradient_function(parameters), function_call_counter

    # These are values deemed good and removed from parameter space for simplification
    MIN_STEP = 1E-8
    BACKTRACK_EPS = 1E-3
    function_call_counter = 0

    if display > 2:
        print("\nValues in bfgs_optimize code:")
        print("\tdimension = %d" % dimension)
        print("\tstep_size = %lg, step_size_adjustment = %lg, reset_when_in_trouble = %s" % (step_size, step_size_adjustment, str(reset_when_in_trouble)))
        print("\trms_gradient = %lg, max_gradient = %lg, max_iterations = %d, MAX_STEP = %lg" % (rms_gradient, max_gradient, max_iterations, MAX_STEP))
        print("\t-------------------")
        print("\treset_step_size = %s, MIN_STEP = %lg, BACKTRACK_EPS = %lg" % (str(reset_step_size), MIN_STEP, BACKTRACK_EPS))
        print("\tfit_rigid = %s\n" % str(fit_rigid))

    # Ensure parameters are in the correct format
    initial_parameters = asarray(initial_parameters).flatten()
    if initial_parameters.ndim == 0:
        initial_parameters.shape = (1,)
    current_parameters = deepcopy(initial_parameters)

    # Initialize inv Hess and Identity matrix
    if start_hess is not None:
        I = start_hess
    else:
        I = eye(len(current_parameters), dtype=float64)
    current_Hessian = I.copy()

    # Get gradient and store your old func_max
    old_fval, current_gradient, function_call_counter = target(current_parameters, function_call_counter)
    grad = current_gradient.flatten().reshape((-1,dimension))
    current_rms_gradient = sqrt((grad**2).sum() / len(grad))
    current_max_gradient = sqrt((grad**2).sum(axis=1).max())

    # Hold original values
    ALPHA_CONST = step_size
    BETA_CONST = step_size_adjustment
    RESET_CONST = reset_step_size

    # Get function to describe linesearch
    if linesearch is 'armijo':
        if display > 1: print("armijo linesearch "),
        def check_backtrack(f1,f0,gk,pk,armijio_line_search_factor,step_size):
            return f1-f0 > armijio_line_search_factor*step_size*dot(gk,pk)
    else:
        if display > 1: print("default linesearch "),
        def check_backtrack(f1,f0,pk,gk,armijio_line_search_factor,step_size):
            return (f1-f0)/(abs(f1)+abs(f0)) > BACKTRACK_EPS

    # Begin Loop
    backtrack, loop_counter, warnflag = 0, 0, 0
    while (current_rms_gradient > rms_gradient) and (function_call_counter < max_iterations) and (current_max_gradient < max_gradient):
        grad = current_gradient.flatten().reshape((-1,dimension))
        current_rms_gradient = sqrt((grad**2).sum() / len(grad))
        current_max_gradient = sqrt((grad**2).sum(axis=1).max())

        # Get your step direction and renorm to remove the effect of current_Hessian not being unit
        step_direction = -dot(current_Hessian, current_gradient).reshape((-1,dimension))
        force_mags = (current_gradient.reshape((-1,dimension))**2).sum(axis=1)
        scalar = sqrt(force_mags / (step_direction**2).sum(axis=1))
        step_direction = (step_direction.T * scalar).T

        # If we are doing unreasonably small step sizes, quit
        step_lengths = sqrt( (step_direction**2).sum(axis=1) ) * step_size
        if max(step_lengths) < MIN_STEP:
            if display > 0: print("Error - Step size unreasonable (%lg)" % abs(max(step_length))),
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
            if display > 0: print("Warning - Setting step to max step size"),
            if reset_when_in_trouble: current_Hessian = I.copy()

        # Hold new parameters
        # We do so because we may not want to keep them if it makes the next step bad.
        new_parameters = current_parameters + step_size * step_direction
        if fit_rigid:
            new_parameters, C, current_Hessian_tmp = align_parameters(new_parameters, [current_gradient, current_parameters], current_Hessian)
            current_gradient_tmp, current_parameters_tmp = C

        # Get the new gradient and check if max has increased
        fval, new_gradient, function_call_counter = target(new_parameters, function_call_counter)

        if target_function is not None and check_backtrack(fval, old_fval, new_gradient, step_direction, armijio_line_search_factor, step_size):
            # Step taken overstepped the minimum.  Lowering step size
            if display > 0:
                print("\tResetting System as %lg > %lg!" % (fval, old_fval))
                print("\talpha: %lg" % step_size),

            step_size *= float64(step_size_adjustment)

            if display > 0: print("-> %lg\n" % step_size)

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
            current_gradient, current_parameters, current_Hessian = current_gradient_tmp, current_parameters_tmp, current_Hessian_tmp

        # Recalculate change_in_parameters to maintain the secant condition
        change_in_parameters = new_parameters - current_parameters
        
        # Store new max value in old_max for future comparison
        if target_function is not None: old_fval = fval

        # Get difference in gradients for further calculations
        change_in_gradient = new_gradient - current_gradient

        try:  # this was handled in numeric, let it remaines for more safety
            rhok = 1.0 / (dot(change_in_gradient, change_in_parameters))
        except ZeroDivisionError:
            rhok = 1000.0
            if display > 1:
                print("Divide-by-zero encountered: rhok assumed large")
        if isinf(rhok):  # this is patch for np
            rhok = 1000.0
            if display > 1:
                print("Divide-by-zero encountered: rhok assumed large")


        # Run BFGS Update for the Inverse Hessian
        A1 = I - change_in_parameters[:, newaxis] * change_in_gradient[newaxis, :] * rhok
        A2 = I - change_in_gradient[:, newaxis] * change_in_parameters[newaxis, :] * rhok
        current_Hessian = dot(A1, dot(current_Hessian, A2)) + (rhok * change_in_parameters[:, newaxis] * change_in_parameters[newaxis, :])

        if display > 1: print("fval %lg" % (fval))

        # Store new parameters, as it has passed the check
        current_parameters = new_parameters
        current_gradient = new_gradient

        # If callback is desired
        if callback is not None:
            callback(current_parameters)

        # Increment the loop counter
        loop_counter += 1

        if (current_rms_gradient <= rms_gradient):
            break
        
    if target_function is not None: fval = old_fval
    else: fval = float('NaN')

    if current_rms_gradient <= gtol: return GTOL_CONVERGENCE
    if loop_counter >= maxiter: return MAXITER_CONVERGENCE
    return FAIL_CONVERGENCE