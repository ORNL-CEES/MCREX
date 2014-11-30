addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'
matrix='laplacian_2d';

addpath(strcat('../utils/model_problems/', matrix))

if strcmp(matrix, 'simple')
    dimen=50;
    A=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
    rhs=ones(dimen,1);
    u=A\rhs;
    
else
     [A, dimen, ~, ~] = mmread('A.mtx');
     rhs=mmread('b.mtx');
     u=mmread('x.mtx');
end

A=sparse(A);

Prec=diag(diag(A));
Prec=sparse(Prec);

H=sparse(speye(size(A))-Prec\A);

%% Numerical setting

numer.eps=10^(-3);
numer.rich_it=300;

%% Statistical setting
stat.nwalks=size(rhs,1);
stat.max_step=1000;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.nchecks=1;
stat.varcut=0.5;
dist=1;

%% Definition of initial and transitional probabilities
[Pb, cdfb, P, cdf]=prob_adjoint(H, rhs, dist);

%% Preconditioning setting

fp.u=u; %reference solution
fp.H=H; %iteration matrix
fp.rhs=rhs;
fp.precond='diag';
%% Monte Carlo Adjoint Method resolution

start=cputime;
%parjob=parpool('local');
%[sol, rel_residual, RES, VAR, DX, NWALKS, tally, iterations]=MCSA_adjoint(fp, dist, P, cdf, numer, stat);
[sol, rel_residual, ~, ~, ~, NWALKS, ~, iterations]=MCSA_adjoint(fp, dist, P, cdf, numer, stat);
%delete(parjob);
finish=cputime;

bar(NWALKS);

save(strcat('../results/MCSA_adjoint/MCSA_adjoint_test_', matrix, '_p=', num2str(dist)))