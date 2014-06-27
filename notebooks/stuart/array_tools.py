##---------------------------------------------------------------------------##
## array_tools.py
## Tools for operating on arrays.
##---------------------------------------------------------------------------##

import numpy

##---------------------------------------------------------------------------##
## Make a Poisson operator with cell centered differences.
##---------------------------------------------------------------------------##
def makePoissonOperator( grid_size, h, D ):
    A = numpy.zeros((grid_size,grid_size))
    A[0][0] = 3.0/(h*h) + D
    A[0][1] = -1.0/(h*h)
    for i in xrange(1,grid_size-1):
        A[i][i-1] = -1.0/(h*h)
        A[i][i] = 2.0/(h*h) + D
        A[i][i+1] = -1.0/(h*h)
    A[grid_size-1][grid_size-2] = -1.0/(h*h)
    A[grid_size-1][grid_size-1] = 3.0/(h*h) + D
    return A

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
## Matrix-matrix multiplication
##---------------------------------------------------------------------------##
def mmMultiply( A, B ):
    num_rows = len(A)
    num_cols = len(B[0])
    N = len(B)
    C = numpy.zeros( (num_rows,num_cols) )
    for i in xrange( num_rows ):
        for j in xrange( num_cols ):
            for k in xrange( N ):
                C[i][j] = C[i][j] + A[i][k]*B[k][j]
    return C

##---------------------------------------------------------------------------##
## Vector update x = x+y
##---------------------------------------------------------------------------##
def updateVector( x, y ):
    size = len(x)
    for i in xrange(size):
        x[i] = x[i] + y[i]
    return x

##---------------------------------------------------------------------------##
## end array_tools.py
##---------------------------------------------------------------------------##

