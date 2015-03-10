addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'; 'sp1_ainv': 'SPN_ainv';

matrix='thermal_eq_diff';

[H,rhs, precond, Prec]=fixed_point(matrix);
%% Numerical setting

numer.eps=10^(-3);
numer.rich_it=300;

%% Statistical setting

stat.nwalks=dimen;
stat.max_step=20;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.varcut=0.1;
dist=1;

%% Definition of initial and transitional probabilities

[Pb, cdfb, P, cdf]=prob_adjoint(H, rhs, dist);

%% Preconditioning setting

fp.u=u; %reference solution
fp.H=H; %iteration matrix
fp.rhs=rhs;
fp.precond='diag';
%% Monte Carlo Adjoint Method resolution

[sol, rel_residual, var, VAR, DX, NWALKS, tally, iterations]=SEQ_adjoint(fp, dist, P, cdf, numer, stat);

conf=0.05;

plot(sol, '*');
hold on
plot(sol-var*norminv(1-conf/2,0,1), 'g*');
plot(sol+var*norminv(1-conf/2,0,1), 'g*');
plot(u, 'r*');


save(strcat('../results/SEQ_adjoint/SEQ_adjoint_test_', matrix, '_p=', num2str(dist)));