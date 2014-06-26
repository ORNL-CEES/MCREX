##---------------------------------------------------------------------------##
## mc_tools.py
## Tools for performing Monte Carlo calculations.
##---------------------------------------------------------------------------##

import math
import random
import mc_data
import numpy
import multilevel_tools

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
    while weight > w_f:
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
## end mc_tools.py
##---------------------------------------------------------------------------##


