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
precond='triblock';

%Sequential Monte Carlo or Monte Carlo Synthetic Acceleration
%possible choices: 'SEQ', 'MCSA'
scheme='MCSA';
method='adjoint';

[fixed_point]=iteration_matrix(precond, D, u, rhs, G);

%spy_matrices(fixed_point);

%% Numerical setting

numer.eps=10^(-7);
numer.rich_it=100000;%maximal number of Richardson iterations

%% Statistical setting

stat.nwalks=100;
stat.max_step=10;
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


%% Post-processing rendering

if reac==0
    if strcmp(scheme, 'SEQ')
        save(strcat('../results/diffusion_problem/Sequential_MC/laplace_adjoint_', scheme, '_p=_', num2str(dist), '_', precond));
    elseif strcmp(scheme, 'MCSA')
        save(strcat('../results/diffusion_problem/MCSA/laplace_adjoint_', scheme, '_p=_', num2str(dist), '_', precond));
    end
else
    if strcmp(scheme, 'SEQ')
        save(strcat('../results/DR_problem/Sequential_MC/laplace_adjoint_', scheme, '_p=_', num2str(dist), '_', precond));
    elseif strcmp(scheme, 'MCSA')
        save(strcat('../results/DR_problem/MCSA/laplace_adjoint_', scheme, '_p=_', num2str(dist), '_', precond));
    end      
end


