addpath('../core')
addpath('../utils')

matrix = 'parabolic_ifiss';

addpath(strcat('../utils/model_problems/', matrix))

if strcmp(matrix, 'parabolic_ifiss')
    load('A.mat');
else
    [A, dimen, ~, ~] = mmread('A.mtx');
end

iperm = load('../utils/metis/graph.txt.iperm');

A = A(iperm+1,iperm+1);

A = tril(A,0);
u=ones(size(A,1),1);
rhs = A*u;
%% 

%algebraic splitting: 'diag', 'gs', 'triblock', 'trisplit', 'alternating'
precond='block_diag';
size_block=1;

assert(mod(size(A,1),size_block)==0, 'the size of the block is not a divisor of the size of the matrix');

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

stat.nwalks=10^4;
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
[sol, rel_residual, VAR, ~, ~, NWALKS, ~, iterations, ~]=system_solver(scheme, method, fixed_point, dist, P, cdf, numer, stat);
%delete(parjob);
finish=cputime;
%% Solver call

start_richardson=cputime;
[sol_richardson, iterations_richardson, relres]=Richardson(fixed_point, numer);
finish_richardson=cputime;

save(strcat('../results/MCSA_adjoint/MCSA_adjoint_test_', matrix, '_p=', num2str(dist)))