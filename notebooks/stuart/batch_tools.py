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

    # Solve the least-squares problem with QR factorization and back
    # out the coefficients
    Q,R = numpy.linalg.qr(ZTZ)
    Q_T = array_tools.matrixTranspose(Q)
    P = numpy.dot(Q_T,ZTb)
    beta = numpy.dot( numpy.linalg.inv(R), P )
    alpha = numpy.dot(W,beta)
    alpha[num_batch-1] = alpha[num_batch-1] + 1

    # Combine the batch solutions using the least-squares coefficients. Also
    # compute a result with a uniform set of coefficients to get the result
    # as if there were no batches.
    x_bmr = numpy.zeros(grid_size)
    for b in xrange(num_batch):
        for i in xrange(grid_size):
            x_bmr[i] = x_bmr[i] + alpha[b]*x_batch[b][i]

    # Return the batch result.
    return x_bmr

##---------------------------------------------------------------------------##
## end batch_tools.py
##---------------------------------------------------------------------------##


