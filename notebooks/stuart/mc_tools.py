##---------------------------------------------------------------------------##
## mc_tools.py
## Tools for performing Monte Carlo calculations.
##---------------------------------------------------------------------------##

import math
import random
import mc_data
import numpy
import multilevel_tools
import array_tools
import batch_tools

##---------------------------------------------------------------------------##
## Given a random number, get an initial state by sampling the source CDF.
##---------------------------------------------------------------------------##
def sampleSourceCDF( rand, source_c ):
    size = len(source_c)
    for i in xrange(size):
        if rand < source_c[i]:
            return i
    return -1

##---------------------------------------------------------------------------##
## Given a state, random number, and matrix CDF, get the new state.
##---------------------------------------------------------------------------##
def sampleMatrixCDF( state, rand, C ):
    row_size = len(C)
    for j in xrange(row_size):
        if rand < C[state][j]:
            return j
    return -1

##---------------------------------------------------------------------------##
## Process a history.
##---------------------------------------------------------------------------##
def doOneHistory( b, source_c, starting_weight, C, W, x, sigma, w_f ):
    size = len(x)
    tally = numpy.zeros(size)
    state = sampleSourceCDF( random.random(), source_c )
    weight = starting_weight*numpy.sign(b[state])
    while abs(weight) > w_f:
        tally[state] = tally[state] + weight
        new_state = sampleMatrixCDF( state, random.random(), C )
        weight = weight * W[state][new_state]
        state = new_state
    for i in xrange(size):
        x[i] = x[i] + tally[i]
        sigma[i] = sigma[i] + tally[i]*tally[i]
    return (x, sigma)

##---------------------------------------------------------------------------##
## Process a multilevel history.
##---------------------------------------------------------------------------##
def doOneMultilevelHistory( b, source_c, starting_weight, C, W, x, sigma, w_f ):
    size = len(x)
    tally_f = numpy.zeros(size)
    state = sampleSourceCDF( random.random(), source_c )
    weight = starting_weight*numpy.sign(b[state])
    while weight > w_f:
        tally_f[state] = tally_f[state] + weight
        new_state = sampleMatrixCDF( state, random.random(), C )
        weight = weight * W[state][new_state]
        state = new_state
    tally_c = multilevel_tools.PRCC( tally_f )
    for i in xrange(size):
        x[i] = x[i] + tally_f[i] - tally_c[i]
        sigma[i] = sigma[i] + (tally_f[i]-tally_c[i])*(tally_f[i]-tally_c[i])
    return (x, sigma)

##---------------------------------------------------------------------------##
## Solve a linear problem using Monte Carlo
##---------------------------------------------------------------------------##
def monteCarloSolve( A, b, w_c, np ):
    H, C, W = mc_data.makeMonteCarloHCW( A )
    source_c, starting_weight = mc_data.makeSourceCDF( b )
    w_f = w_c * starting_weight
    x = numpy.zeros( len(b) )
    sigma = numpy.zeros( len(b) )
    for i in xrange(np):
        x, sigma = doOneHistory( b, source_c, starting_weight, C, W, x, sigma, w_f )
    for i in xrange(len(x)):
        x[i] = x[i] / np
        sigma[i] = (sigma[i] / np - x[i]*x[i])/np
    return (x, sigma)

##---------------------------------------------------------------------------##
## Solve a linear problem using a multilevel estimator
##---------------------------------------------------------------------------##
def multilevelMonteCarloSolve( A, b, w_c, np ):
    H, C, W = mc_data.makeMonteCarloHCW( A )
    source_c, starting_weight = mc_data.makeSourceCDF( b )
    w_f = w_c * starting_weight
    x = numpy.zeros( len(b) )
    sigma = numpy.zeros( len(b) )
    for i in xrange(np):
        x, sigma = doOneMultilevelHistory( b, source_c, starting_weight, C, W, x, sigma, w_f )
    for i in xrange(len(x)):
        x[i] = x[i] / np
        sigma[i] = (sigma[i] / np - x[i]*x[i])/np
    return (x, sigma)

##---------------------------------------------------------------------------##
## Solve a linear problem using Monte Carlo and return the result in batches.
##---------------------------------------------------------------------------##
def batchMonteCarloSolve( A, b, w_c, np ):
    H, C, W = mc_data.makeMonteCarloHCW( A )
    source_c, starting_weight = mc_data.makeSourceCDF( b )
    w_f = w_c * starting_weight
    x = []
    tally = numpy.zeros( len(b) )
    sigma = numpy.zeros( len(b) )
    for i in xrange(np):
        dummy = numpy.zeros( len(b) )
        tally, sigma = doOneHistory( b, source_c, starting_weight, C, W, dummy, sigma, w_f )
        x.append(tally)
    return (x, sigma)

##---------------------------------------------------------------------------##
## Solve a linear problem with MCSA
##---------------------------------------------------------------------------##
def solveMCSA( A, x, b, tol, max_iter, w_c, np ):
    r = array_tools.computeResidual(A,x,b)
    r_norm = numpy.linalg.norm(r,2)
    b_norm = numpy.linalg.norm(b,2)
    iter = 0
    while ( r_norm/b_norm > tol ) and ( iter < max_iter ):
        x = array_tools.updateVector(x,r)
        r = array_tools.computeResidual(A,x,b)
        delta, sigma = monteCarloSolve( A, r, w_c, np )
        x = array_tools.updateVector(x,delta)
        r = array_tools.computeResidual(A,x,b)
        r_norm = numpy.linalg.norm(r,2)
        iter = iter + 1
        print iter, ":", r_norm / b_norm
    return x

##---------------------------------------------------------------------------##
## Solve a linear problem with MCSA and batch minimization
##---------------------------------------------------------------------------##
def solveMRBMCSA( A, x, b, tol, max_iter, w_c, np, num_batch ):
    r = array_tools.computeResidual(A,x,b)
    r_norm = numpy.linalg.norm(r,2)
    b_norm = numpy.linalg.norm(b,2)
    iter = 0
    while ( r_norm/b_norm > tol ) and ( iter < max_iter ):
        x = array_tools.updateVector(x,r)
        r = array_tools.computeResidual(A,x,b)
        delta_sample, sigma = batchMonteCarloSolve( A, r, w_c, np )
        delta = batch_tools.computeMRB( A, delta_sample, r, num_batch)
        x = array_tools.updateVector(x,delta)
        r = array_tools.computeResidual(A,x,b)
        r_norm = numpy.linalg.norm(r,2)
        iter = iter + 1
        print iter, ":", r_norm / b_norm
    return x

##---------------------------------------------------------------------------##
## Process a history for a polynomial with the expected value estimator.
##---------------------------------------------------------------------------##
def doOnePolyHistory( A, b, source_c, starting_weight, C, W, x, rank ):
    size = len(x[0])
    state = sampleSourceCDF( random.random(), source_c )
    weight = starting_weight*numpy.sign(b[state])
    for q in xrange(1,rank):
        for j in xrange(size):
            x[q][j] = x[q][j] + A[j][state]*weight
        new_state = sampleMatrixCDF( state, random.random(), C )
        weight = weight * W[state][new_state]
        state = new_state
    return x

##---------------------------------------------------------------------------##
## Stochastically construct a matrix polynomial b + Ab + (A^2)b + ... 
##---------------------------------------------------------------------------##
def monteCarloMatrixPolynomial( A, b, np, rank ):
    P = mc_data.makeProbabilityMatrix( A )
    C = mc_data.makeCDFMatrix( P )
    W = mc_data.makeWeightMatrix( A, P )
    source_c, starting_weight = mc_data.makeSourceCDF( b )
    x = numpy.zeros( (rank,len(b)) )
    x[0] = b
    for i in xrange(np):
        x = doOnePolyHistory( A, b, source_c, starting_weight, C, W, x, rank )
    for q in xrange(rank):
        for i in xrange(len(x[0])):
            x[q][i] = x[q][i] / np
    return x

##---------------------------------------------------------------------------##
## end mc_tools.py
##---------------------------------------------------------------------------##


