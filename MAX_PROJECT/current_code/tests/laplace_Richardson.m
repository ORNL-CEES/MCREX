addpath('../core')
addpath('../utils')

%parallel
%poolobj=parpool('local');

shape = 'S'; % Other possible shapes include L,S,N,C,D,A,H,B
n=100;

% creation of te grid
G=gridgen(shape, n);

reac=0.1;
% computation of the solution with \
[u, D, rhs]=laplace(shape, G, reac);

%algebraic splitting: 'diag', 'gs', 'triblock', 'trisplit', 'alternating'
precond='diag';

%Sequential Monte Carlo or Monte Carlo Synthetic Acceleration
%possible choices: 'SEQ', 'MCSA'
scheme='MCSA';
method='adjoint';

[fixed_point]=iteration_matrix(precond, D, u, rhs, G);

%spy_matrices(fixed_point);

%% Numerical setting

numer.eps=10^(-7);
numer.rich_it=100000;%maximal number of Richardson iterations
%% Solver call

start=cputime;
[sol, iterations, relres]=Richardson(fixed_point, numer);
finish=cputime;
