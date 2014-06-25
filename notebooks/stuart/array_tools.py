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
## Get the inverse diagonal of a matrix.
##---------------------------------------------------------------------------##
def getInvDiag( A ):
    size = len(A)
    diag = numpy.zeros(size)
    for i in xrange(size):
        diag[i] = 1.0 / A[i][i]
    return diag

##---------------------------------------------------------------------------##
## Left scale a matrix by a vector.
##---------------------------------------------------------------------------##
def leftScaleMatrix( A, v ):
    size = len(A)
    for i in xrange(size):
        for j in xrange(size):
            A[i][j] = A[i][j] * v[i]
    return A

##---------------------------------------------------------------------------##
## Scale a vector by a vector.
##---------------------------------------------------------------------------##
def scaleVector( x, v ):
    size = len(x)
    for i in xrange(size):
        x[i] = x[i] * v[i]
    return x

##---------------------------------------------------------------------------##
## end array_tools.py
##---------------------------------------------------------------------------##

