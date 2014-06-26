##---------------------------------------------------------------------------##
## multilevel_tools.py
## Tools for multilevel Monte Carlo calculations.
##---------------------------------------------------------------------------##

import array_tools
import numpy

##---------------------------------------------------------------------------##
## Build a prolongation operator to a given matrix.
##---------------------------------------------------------------------------##
def buildP( A ):
    num_rows = len(A)
    num_cols = num_rows / 2
    P = numpy.zeros( (num_rows,num_cols) )
    for i in xrange( num_cols ):
        P[2*i][i] = 1.0
        P[2*i+1][i] = 1.0
    return P

##---------------------------------------------------------------------------##
## Build a restriction operator from a given matrix.
##---------------------------------------------------------------------------##
def buildR( A ):
    num_cols = len(A)
    num_rows = num_cols / 2
    R = numpy.zeros( (num_rows,num_cols) )
    for j in xrange( num_rows ):
        R[j][2*j] = 0.5
        R[j][2*j+1] = 0.5
    return R

##---------------------------------------------------------------------------##
## Build a coarse level operator from a fine level operator.
##---------------------------------------------------------------------------##
def buildRAP( A ):
    P = buildP( A )
    R = buildR( A )
    AP = array_tools.mmMultiply( A, P )
    RAP = array_tools.mmMultiply( R, AP )
    return RAP

##---------------------------------------------------------------------------##
## Prolongate a result from a cell-centered coarse grid to a
## cell-centered fine grid.
## ---------------------------------------------------------------------------##
def prolongateCC( u_c ):
    coarse_size = len(u_c)
    fine_size = 2*coarse_size
    u_f = numpy.zeros( fine_size )
    for i in xrange(coarse_size):
        u_f[2*i] = u_c[i]
        u_f[2*i+1] = u_c[i]
    return u_f

##---------------------------------------------------------------------------##
## Restrict a result from a cell-centered fine grid to a cell-centered
## coarse grid.
## ---------------------------------------------------------------------------##
def restrictCC( u_f ):
    fine_size = len(u_f)
    coarse_size = fine_size / 2
    u_c = numpy.zeros( coarse_size )
    for i in xrange(coarse_size):
        u_c[i] = 0.5*(u_f[2*i] + u_f[2*i+1])
    return u_c

##---------------------------------------------------------------------------##
## Restrict to a coarse grid and then immediately prolongate back on a
## cell-centered grid.
## ---------------------------------------------------------------------------##
def PRCC( u ):
    u_c = restrictCC( u )
    u_f = prolongateCC( u_c )
    return u_f

##---------------------------------------------------------------------------##
## end multilevel_tools.py
##---------------------------------------------------------------------------##

