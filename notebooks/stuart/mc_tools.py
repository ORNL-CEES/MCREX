##---------------------------------------------------------------------------##
## mc_tools.py
## Tools for performing Monte Carlo calculations.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Given a random number, get an initial state by sampling the source CDF.
##---------------------------------------------------------------------------##
def sampleSourceCDF( rand, source_c ):
    size = len(source_c)
    for i in xrange(size):
        if rand >= source_c[i]:
            return i
    return -1

##---------------------------------------------------------------------------##
## Given a state, random number, and matrix CDF, get the new state.
##---------------------------------------------------------------------------##
def sampleMatrixCDF( in_state, rand, C ):
    row_size = len(C)
    for j in xrange(row_size):
        if rand >= C[in_state][j]:
            return j
    return -1

##---------------------------------------------------------------------------##
## end mc_tools.py
##---------------------------------------------------------------------------##


