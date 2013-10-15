##---------------------------------------------------------------------------##
## matrix.py
##---------------------------------------------------------------------------##

import os, sys

##---------------------------------------------------------------------------##
## Build an 7-point stencil matrix
##---------------------------------------------------------------------------##

Nx = 10
Ny = 10
Nz = 5
N  = Nx*Ny*Nz

A = []

for k in xrange(Nz):
    for j in xrange(Ny):
        for i in xrange(Nx):
            row = [0] * N
            
            px = i + 1 + Nx * (j + k * Ny)
            mx = i - 1 + Nx * (j + k * Ny)
            py = i + Nx * (j + 1 + k * Ny)
            my = i + Nx * (j - 1 + k * Ny)
            pz = i + Nx * (j + (k + 1) * Ny)
            mz = i + Nx * (j + (k - 1) * Ny)

            c  = i + Nx * (j + k * Ny)
            
            row[c] = 1

            if (i > 0) : row[mx] = 1
            if (j > 0) : row[my] = 1
            if (k > 0) : row[mz] = 1

            if (i < Nx - 1) : row[px] = 1
            if (j < Ny - 1) : row[py] = 1
            if (k < Nz - 1) : row[pz] = 1

            A.append(row)
            
for j in xrange(N):
    for i in xrange(N):
        print A[i][j],
    print


##---------------------------------------------------------------------------##
## end of matrix.py
##---------------------------------------------------------------------------##


