addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'; 'sp1_ainv': 'SPN_ainv';

matrix='SPN_ainv';
[H,rhs, precond, Prec]=fixed_point(matrix);

%% Numerical setting

numer.eps=10^(-6);
numer.rich_it=3000;

%% Statistical setting

stat.nwalks=100;
stat.max_step=10;
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
fp.precond=precond;

%% Monte Carlo Adjoint Method resolution

start=cputime;
%parjob=parpool('local');
% [sol, rel_residual, RES, VAR, DX, NWALKS, tally, iterations]=MCSA_adjoint(fp, dist, P, cdf, numer, stat);
[sol, rel_residual, ~, ~, ~, NWALKS, ~, iterations]=MCSA_adjoint(fp, dist, P, cdf, numer, stat);
%delete(parjob);
finish=cputime;

if strcmp(matrix, 'sp1_shift') || strcmp(matrix, 'sp3_shift') || strcmp(matrix, 'SPN_shift') || strcmp(matrix, 'sp5_shift') ...
        || strcmp(matrix, 'sp1_ainv') || strcmp(matrix, 'SPN_ainv')
    sol=Prec\sol;
end

bar(NWALKS);

save(strcat('../results/MCSA_adjoint/MCSA_adjoint_test_', matrix, '_p=', num2str(dist)))