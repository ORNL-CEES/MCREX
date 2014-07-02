##---------------------------------------------------------------------------##
## variance_reduction.py
## Variance reduction for monte carlo calculations.
##---------------------------------------------------------------------------##

import math
import random
import numpy

##---------------------------------------------------------------------------##
## For a given importance map, window size, and survival weight
## constant, compute the upper and lower bounds of the weight windows
## and the survival weights.
## ---------------------------------------------------------------------------##
def computeWeightWindows( I,  c_u, c_s ):
    R = numpy.linalg.norm( I, 1 )
    size = len(I)
    w_l = numpy.zeros( size )
    w_u = numpy.zeros( size )
    w_s = numpy.zeros( size )
    for i in xrange(size):
        w_l[i] = 2.0 * R / ( abs(I[i])*(c_u + 1 ) )
        w_u[i] = c_u * w_l[i]
        w_s[i] = c_s * w_l[i]
    return (w_l,w_u,w_s)

##---------------------------------------------------------------------------##
## For a given sample weight and survivial weight, determine if
## roulette should be performed. Return true if we should kill the
## history.
## ---------------------------------------------------------------------------##
def roullette( w, w_s ):
    rn = random.random()
    return rn > (w/w_s)

##---------------------------------------------------------------------------##
## For a given sample, split it into a number of samples with weight
## in the center of the weight window. Return the number of new
## samples and their weight.
## ---------------------------------------------------------------------------##
def split( w, w_l, w_u ):
    w_mid = (w_u - w_l) / 2.0
    n_split = math.floor( w / w_mid )
    w_split = w / n_split
    if ( w_split < w_l or w_split > w_u ):
        print "Split weight outside of bounds!"
        raise
    return n_split, w_split

##---------------------------------------------------------------------------##
## end variance_reduction.py
##---------------------------------------------------------------------------##
