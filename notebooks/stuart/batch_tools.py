##---------------------------------------------------------------------------##
## batch_tools.py
## Tools for batch Monte Carlo.
##---------------------------------------------------------------------------##

import numpy
import array_tools

##---------------------------------------------------------------------------##
## Given a set of individual sample results, compute the
## minimized-residual batch result for the given batch size.
## ---------------------------------------------------------------------------##
def computeMRB( A, x_sample, f, num_batch ):

    # Combine individual sample results into batches.
    grid_size = len(x_sample[0])
    num_sample = len(x_sample)
    sample_per_batch = num_sample / num_batch
    x_batch = numpy.zeros( (num_batch,grid_size) )
    batch_norm = 1.0 / sample_per_batch
    for i in xrange(num_batch):
        for j in xrange(grid_size):
            for k in xrange(sample_per_batch):
                batch_idx = sample_per_batch*i + k
                x_batch[i][j] = x_batch[i][j] + batch_norm*x_sample[batch_idx][j]

    # If there was a single batch return the result.
    if ( 1 == num_batch ):
        return x_batch[0]

    # Make the elimination matrix, W, and its transpose
    W = numpy.zeros( (num_batch,num_batch-1) )
    W[0][0] = -1
    for i in xrange(1,num_batch-1):
        W[i][i-1] = 1
        W[i][i] = -1
    W[num_batch-1][num_batch-2] = 1
    W_T = array_tools.matrixTranspose(W)

    # Make the least-squares problem operator.
    V_T = []
    for i in xrange(num_batch):
        Ax_b = numpy.dot(A,x_batch[i])
        V_T.append( Ax_b )
    Z_T = numpy.dot(W_T,V_T)
    V = array_tools.matrixTranspose(V_T)
    Z = numpy.dot(V,W)
    ZTZ = numpy.dot(Z_T,Z)

    # Make the least-squares problem RHS
    b = numpy.zeros(grid_size)
    for i in xrange(grid_size):
        b[i] = f[i] - V[i][num_batch-1]
    ZTb = numpy.dot(Z_T,b)

    # Solve the least-squares problem with QR decomposition and back
    # out the coefficients
    Q,R = numpy.linalg.qr(ZTZ)
    Q_T = array_tools.matrixTranspose(Q)
    P = numpy.dot(Q_T,ZTb)
    beta = numpy.dot( numpy.linalg.inv(R), P )
    alpha = numpy.dot(W,beta)
    alpha[num_batch-1] = alpha[num_batch-1] + 1

    # Combine the batch solutions using the least-squares coefficients.
    x_bmr = numpy.zeros(grid_size)
    for b in xrange(num_batch):
        for i in xrange(grid_size):
            x_bmr[i] = x_bmr[i] + alpha[b]*x_batch[b][i]

    # Return the batch result.
    return x_bmr

##---------------------------------------------------------------------------##
## Given a set of basis vectors from which to extract a correction,
## find their linear combination that minimizes the linear problem
## residual by solving an linearly constrained linear least squares
## problem.
## ---------------------------------------------------------------------------##
def computeConstrainedLeastSquares( A, basis, f, x_0 ):

    num_basis = len(basis)
    grid_size = len(basis[0])

    # If there was a single basis vector return the result.
    if ( 1 == num_basis ):
        return basis[0]

    # Make the elimination matrix, W, and its transpose
    W = numpy.zeros( (num_basis,num_basis-1) )
    W[0][0] = -1
    for i in xrange(1,num_basis-1):
        W[i][i-1] = 1
        W[i][i] = -1
    W[num_basis-1][num_basis-2] = 1
    W_T = array_tools.matrixTranspose(W)

    # Make the least-squares problem operator.
    V_T = []
    for i in xrange(num_basis):
        Ax_b = numpy.dot(A,basis[i])
        V_T.append( Ax_b )
    V = array_tools.matrixTranspose(V_T)
    Z_T = numpy.dot(W_T,V_T)
    Z = array_tools.matrixTranspose(Z_T)
    ZTZ = numpy.dot(Z_T,Z)

    # Make the least-squares problem RHS
    b = numpy.zeros(grid_size)
    Ax_0 = numpy.dot(A,x_0)
    for i in xrange(grid_size):
        b[i] = f[i] - V[i][num_basis-1] - Ax_0[i]
    ZTb = numpy.dot(Z_T,b)

    # Solve the least-squares problem with QR decomposition and back
    # out the coefficients
    Q,R = numpy.linalg.qr(ZTZ)
    Q_T = array_tools.matrixTranspose(Q)
    P = numpy.dot(Q_T,ZTb)
    beta = numpy.dot( numpy.linalg.inv(R), P )
    alpha = numpy.dot(W,beta)
    alpha[num_basis-1] = alpha[num_basis-1] + 1

    # Combine the basis vectors using the least-squares coefficients.
    x = numpy.zeros(grid_size)
    for b in xrange(num_basis):
        for i in xrange(grid_size):
            x[i] = x[i] + alpha[b]*basis[b][i]

    # Return the combined result.
    return x

##---------------------------------------------------------------------------##
## Given a set of basis vectors, find their linear combination that
## minimizes the linear problem residual by solving an unconstrained
## linear least squares problem.
## ---------------------------------------------------------------------------##
def computeLeastSquares( A, basis, f, x_0 ):

    num_basis = len(basis)
    grid_size = len(basis[0])

    # If there was a single basis vector return the result.
    if ( 1 == num_basis ):
        return basis[0]

    # Make the least-squares problem operator.
    V_T = []
    for i in xrange(num_basis):
        Ax_b = numpy.dot(A,basis[i])
        V_T.append( Ax_b )
    V = array_tools.matrixTranspose(V_T)
    VTV = numpy.dot(V_T,V)

    # Make the least-squares problem RHS
    b = numpy.zeros(grid_size)
    Ax_0 = numpy.dot(A,x_0)
    for i in xrange(grid_size):
        b[i] = f[i] - V[i][num_basis-1] - Ax_0[i]
    VTb = numpy.dot(V_T,b)

    # Solve the least-squares problem with QR decomposition and back
    # out the coefficients
    Q,R = numpy.linalg.qr(VTV)
    Q_T = array_tools.matrixTranspose(Q)
    P = numpy.dot(Q_T,VTb)
    alpha = numpy.dot( numpy.linalg.inv(R), P )

    # Combine the basis vectors using the least-squares coefficients.
    x = numpy.zeros(grid_size)
    for b in xrange(num_basis):
        for i in xrange(grid_size):
            x[i] = x[i] + alpha[b]*basis[b][i]

    # Return the combined result.
    return x

##---------------------------------------------------------------------------##
## end batch_tools.py
##---------------------------------------------------------------------------##


