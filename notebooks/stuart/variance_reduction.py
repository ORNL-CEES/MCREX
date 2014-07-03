##---------------------------------------------------------------------------##
## variance_reduction.py
## Variance reduction for monte carlo calculations.
##---------------------------------------------------------------------------##

import math
import random
import numpy
import mc_data
import mc_tools
import array_tools

##---------------------------------------------------------------------------##
## For a given importance map, source function, window size, and
## survival weight constant, compute the upper and lower bounds of the
## weight windows and the survival weights.
## ---------------------------------------------------------------------------##
def computeWeightWindows( I, R, c_u ):
    size = len(I)
    w_l = numpy.zeros( size )
    w_u = numpy.zeros( size )
    w_s = numpy.zeros( size )
    for i in xrange(size):
        w_l[i] = 2.0 * R / (abs(I[i])*(c_u + 1 ) )
        w_u[i] = c_u * w_l[i]
        w_s[i] = (w_l[i] + w_u[i]) / 2.0
    return (w_l,w_u,w_s)

##---------------------------------------------------------------------------##
## For a given sample weight and survivial weight, determine if
## roulette should be performed. Return true if we should kill the
## history.
## ---------------------------------------------------------------------------##
def roulette( w, w_s ):
    rn = random.random()
    return rn > (w/w_s)

##---------------------------------------------------------------------------##
## For a given sample, split it into a number of samples with weight
## in the center of the weight window. Return the number of new
## samples and their weight.
## ---------------------------------------------------------------------------##
def split( w, w_l, w_u, w_s ):
    n_split = int(math.floor( w / w_s ))
    w_split = w / n_split
    if ( w_split < w_l or w_split > w_u ):
        print "ERROR: Split weight outside of bounds!"
    return (n_split, w_split)

##---------------------------------------------------------------------------##
## Do one history with the adjoint method and weight windows. This
## includes all histories created from splitting.
## ---------------------------------------------------------------------------##
def doOneHistory( b, source_c, starting_weight, C, W, x, sigma, w_f, w_l, w_u, w_s ):
    history_stack = []
    weight_stack = []
    hw_stack = []
    size = len( x )
    tally = numpy.zeros(size)
    state = mc_tools.sampleSourceCDF( random.random(), source_c )
    history_stack.append(state)
    weight_stack.append( starting_weight*numpy.sign(b[state]) )
    hw_stack.append(starting_weight)

    # Run all histories in the stack.
    while len(history_stack) > 0:
        state = history_stack.pop()
        weight = weight_stack.pop()
        hw = hw_stack.pop()

        # Run histories until the cutoff.
        while abs(weight) > w_f:

            # Process a transition
            tally[state] = tally[state] + (hw/starting_weight)*weight
            new_state = mc_tools.sampleMatrixCDF( state, random.random(), C )
            weight = weight * W[state][new_state]
            state = new_state

            # Roulette if below weight window.
            if (hw < w_l[state]):

                # Set the weight to zero to kill the history.
                if ( roulette(hw,w_s[state]) ):
                    weight = 0.0

                # If it survived set the weight to the survival weight.
                else:
                    hw = w_s[state]

            # Split if above weight window
            elif (hw > w_u[state]):

                # Compute the number of new histories and push them
                # onto the stack.
                n_split, w_split = split( hw, w_l[state], w_u[state] )
                for j in xrange(n_split):
                    history_stack.append(state)
                    weight_stack.append(weight)
                    hw_stack.append(w_split)

                # End the current history in leiu of the new histories
                # from the split.
                weight = 0.0
        
    # Combine the tally results for the mean and variance.
    for i in xrange(size):
        x[i] = x[i] + tally[i]
        sigma[i] = sigma[i] + tally[i]*tally[i]

    # Return the mean and variance estimates.
    return (x, sigma)

##---------------------------------------------------------------------------##
## Solve a linear problem using nonanalog Monte Carlo
##---------------------------------------------------------------------------##
def monteCarloSolve( A, b, w_c, np, I, R, c_u ):
    H, C, W = mc_data.makeMonteCarloHCW( A )
    source_c, starting_weight = mc_data.makeSourceCDF( b )
    w_f = w_c * starting_weight
    w_l,w_u,w_s = computeWeightWindows( I, R, c_u )
    x = numpy.zeros( len(b) )
    sigma = numpy.zeros( len(b) )

    for i in xrange(np):
        x, sigma = doOneHistory( b, source_c, starting_weight, C, W, x, sigma, w_f, w_l, w_u, w_s )
    for i in xrange(len(x)):
        x[i] = x[i] / np
        sigma[i] = (sigma[i] / np - x[i]*x[i])/np
    return (x, sigma)

##---------------------------------------------------------------------------##
## FW-Cadis solve
##---------------------------------------------------------------------------##
def fwcadis( A, b, w_c, np, c_u ):
    x_f = numpy.linalg.solve(A,b)
    A_T = array_tools.matrixTranspose( A )
    d = numpy.zeros(len(x_f))
    for i in xrange(len(x_f)):
        d[i] = 1.0 / x_f[i]
    y = numpy.linalg.solve(A_T,d)
    R = numpy.dot(y,b)
    print y
    x, sigma = monteCarloSolve( A, b, w_c, np, y, R, c_u )
    return (x,sigma)

##---------------------------------------------------------------------------##
## Solve a linear problem with MCSA
##---------------------------------------------------------------------------##
def solveMCSA( A, x, b, tol, max_iter, w_c, np, c_u ):
    r = array_tools.computeResidual(A,x,b)
    r_norm = numpy.linalg.norm(r,2)
    b_norm = numpy.linalg.norm(b,2)
    iter = 0
    while ( r_norm/b_norm > tol ) and ( iter < max_iter ):
        x = array_tools.updateVector(x,r)
        r = array_tools.computeResidual(A,x,b)
        response = numpy.linalg.norm(r,1)
        delta, sigma = monteCarloSolve( A, r, w_c, np, r, response, c_u )
        x = array_tools.updateVector(x,delta)
        r = array_tools.computeResidual(A,x,b)
        r_norm = numpy.linalg.norm(r,2)
        iter = iter + 1
        print iter, ":", r_norm / b_norm
    return x

##---------------------------------------------------------------------------##
## end variance_reduction.py
##---------------------------------------------------------------------------##
