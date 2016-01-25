addpath('../core')
addpath('../utils')

matrix = 'laplace_nilpotent';


shape = 'S'; % Other possible shapes include L,S,N,C,D,A,H,B
n=32;

% creation of the grid
G=gridgen(shape, n);

reac = 0;
reordering = 1;

% computation of the solution with \
[~, A, ~]=laplace(shape, G, reac);

if reordering
    ci = scomponents(A);
    [B,I] = sort(ci, 'descend');
    A=A(I,I);
end
    

A = tril(A,0);
u = ones(size(A,1),1);

if amd_reordering
    rhs = D * A * u;   
else
    rhs = A * u;
end
%% 

%algebraic splitting: 'diag', 'gs', 'triblock', 'trisplit', 'alternating'
precond='no';
size_block = 1;

%Sequential Monte Carlo or Monte Carlo Synthetic Acceleration
%possible choices: 'SEQ', 'MCSA'
scheme='MCSA';
method='adjoint';

[fixed_point]=iteration_matrix2(precond, A, u, rhs, size_block);

%spy_matrices(fixed_point);

%% Numerical setting

numer.eps=10^(-7);
numer.rich_it=100000;%maximal number of Richardson iterations

%% Statistical setting

stat.nwalks=1000;
stat.max_step=1000;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.nchecks=1;
stat.varcut=0.1;
dist=1;

%% Definition of initial and transitional prbabilities

if ~ strcmp(precond, 'alternating')
    [Pb, cdfb, P, cdf]=prob_adjoint(fixed_point.H, fixed_point.rhs, dist);
    
else
    [Pb.Pb1, cdfb.cdfb1, P.P1, cdf.cdf1]=prob_adjoint(fixed_point.H1, fixed_point.rhs1, dist);
    [Pb.Pb2, cdfb.cdfb2, P.P2, cdf.cdf2]=prob_adjoint(fixed_point.H2, fixed_point.rhs2, dist);
end

%% Solver call

start=cputime;
%parjob=parpool('local');
[sol, rel_residual, ~, ~, ~, NWALKS, ~, iterations, ~]=system_solver(scheme, method, fixed_point, dist, P, cdf, numer, stat);
%delete(parjob);
finish=cputime;
%% Solver call

start_richardson=cputime;
[sol_richardson, iterations_richardson, relres]=Richardson(fixed_point, numer);
finish_richardson=cputime;

save(strcat('../results/MCSA_adjoint/MCSA_adjoint_test_', matrix, '_p=', num2str(dist)))