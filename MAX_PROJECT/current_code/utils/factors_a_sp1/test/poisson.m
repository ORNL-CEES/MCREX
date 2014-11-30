function A = poisson(n)
%POISSON   Block tridiagonal matrix from Poisson's equation (sparse).
%          POISSON(N) is the block tridiagonal matrix of order N^2
%          resulting from discretizing Poisson's equation with the
%          5-point operator on an N-by-N mesh.

%          Reference:
%          G.H. Golub and C.F. Van Loan, Matrix Computations, second edition,
%          Johns Hopkins University Press, Baltimore, Maryland, 1989
%          (Section 4.5.4).

S = tridiag(n,-1,2,-1);
I = speye(n);
A = kron(I,S) + kron(S,I);
function T = tridiag(n, x, y, z)
%TRIDIAG  Tridiagonal matrix (sparse).
%         TRIDIAG(X, Y, Z) is the tridiagonal matrix with subdiagonal X,
%         diagonal Y, and superdiagonal Z.
%         X and Z must be vectors of dimension one less than Y.
%         Alternatively TRIDIAG(N, C, D, E), where C, D, and E are all
%         scalars, yields the Toeplitz tridiagonal matrix of order N
%         with subdiagonal elements C, diagonal elements D, and superdiagonal
%         elements E.   This matrix has eigenvalues (Todd 1977)
%                  D + 2*SQRT(C*E)*COS(k*PI/(N+1)), k=1:N.
%         TRIDIAG(N) is the same as TRIDIAG(N,-1,2,-1), which is
%         a symmetric positive definite M-matrix (the negative of the
%         second difference matrix).

%         References:
%         J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
%           Birkhauser, Basel, and Academic Press, New York, 1977, p. 155.
%         D.E. Rutherford, Some continuant determinants arising in physics and
%           chemistry---II, Proc. Royal Soc. Edin., 63, A (1952), pp. 232-241.

if nargin == 1, x = -1; y = 2; z = -1; end
if nargin == 3, z = y; y = x; x = n; end

x = x(:); y = y(:); z = z(:);   % Force column vectors.

if max( [ size(x) size(y) size(z) ] ) == 1
   x = x*ones(n-1,1);
   z = z*ones(n-1,1);
   y = y*ones(n,1);
else
   [nx, m] = size(x);
   [ny, m] = size(y);
   [nz, m] = size(z);
   if (ny - nx - 1) | (ny - nz -1)
      error('Dimensions of vector arguments are incorrect.')
   end
end

% T = diag(x, -1) + diag(y) + diag(z, 1);  % For non-sparse matrix.
n = max(size(y));
T = spdiags([ [x;0] y [0;z] ], -1:1, n, n);
