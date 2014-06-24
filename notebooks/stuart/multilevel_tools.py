##---------------------------------------------------------------------------##
## multilevel_tools.py
## Tools for multilevel Monte Carlo calculations.
##---------------------------------------------------------------------------##

import numpy

##---------------------------------------------------------------------------##
## Prolongate a result from a cell-centered coarse grid to a
## cell-centered fine grid.
## ---------------------------------------------------------------------------##
def prolongateCC( u_c ):
    coarse_size = len(u_c)
    fine_size = 2*coarse_size
    u_f = numpy.zeros( fine_size )
    for i in xrange(coarse_size):
        u_f[2*i] = 0.5*u_c[i]
        u_f[2*i+1] = 0.5*u_c[i]
    return u_f

##---------------------------------------------------------------------------##
## Restrict a result from a cell-centered fine grid to a cell-centered
## coarse grid.
## ---------------------------------------------------------------------------##
def restrictCC( u_f ):
    fine_size = len(u_f)
    coarse_size = fine_size / 2
    u_c = numpy.zeros( fine_size )
    for i in xrange(coarse_size):
        u_c[i] = u_f[2*i] + u_f[2*i+1]
    return u_c

##---------------------------------------------------------------------------##
## end multilevel_tools.py
##---------------------------------------------------------------------------##

