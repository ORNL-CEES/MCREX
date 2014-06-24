##---------------------------------------------------------------------------##
## mc_tools.py
## Tools for performing Monte Carlo calculations.
##---------------------------------------------------------------------------##

import random
import mc_data
import numpy

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
def doOneHistory( source_c, starting_weight, C, W, x, w_f ):
    state = sampleSourceCDF( random.random(), source_c )
    weight = starting_weight
    while weight > w_f:
        x[state] = x[state] + weight
        new_state = sampleMatrixCDF( state, random.random(), C )
        weight = weight * W[state][new_state]
        state = new_state
    return x

##---------------------------------------------------------------------------##
## Solve a linear problem using Monte Carlo
##---------------------------------------------------------------------------##
def monteCarloSolve( A, b, w_c, np ):
    H, C, W = mc_data.makeMonteCarloHCW( A )
    source_c, starting_weight = mc_data.makeSourceCDF( b )
    w_f = w_c * starting_weight
    x = numpy.zeros( len(b) )
    for i in xrange(np):
        x = doOneHistory( source_c, starting_weight, C, W, x, w_f )
    for i in xrange(len(x)):
        x[i] = x[i] / np
    return x

##---------------------------------------------------------------------------##
## end mc_tools.py
##---------------------------------------------------------------------------##


