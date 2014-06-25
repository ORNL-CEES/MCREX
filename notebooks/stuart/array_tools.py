##---------------------------------------------------------------------------##
## array_tools.py
## Tools for operating on arrays.
##---------------------------------------------------------------------------##

import numpy

##---------------------------------------------------------------------------##
## Compute the relative difference between two arrays.
##---------------------------------------------------------------------------##
def relativeDiff( x, x_ref ):
    size = len(x)
    diff = numpy.zeros(size)
    for i in xrange(size):
        diff[i] = abs( (x[i]-x_ref[i])/x_ref[i] )
    return diff

##---------------------------------------------------------------------------##
## Compute the 2-norm of the relative difference between two arrays.
##---------------------------------------------------------------------------##
def diff2Norm( x, x_ref ):
    diff = relativeDiff( x, x_ref )
    return numpy.linalg.norm( diff, 2 )

##---------------------------------------------------------------------------##
## Compute the inf-norm of the relative difference between two arrays.
##---------------------------------------------------------------------------##
def diffInfNorm( x, x_ref ):
    diff = relativeDiff( x, x_ref )
    return numpy.linalg.norm( diff, inf )

##---------------------------------------------------------------------------##
## Compute the residual of a linear system.
##---------------------------------------------------------------------------##
def computeResidual( A, x, b ):
    size = len(x)
    r = numpy.zeros(size)
    Ax = numpy.dot(A,x)
    for i in xrange(size):
        r[i] = b[i] - Ax[i]
    return r

##---------------------------------------------------------------------------##
## end array_tools.py
##---------------------------------------------------------------------------##

