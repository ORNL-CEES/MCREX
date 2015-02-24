addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'; 'sp1_ainv': 'SPN_ainv';

matrix='SPN_ainv';
[H,rhs, precond, Prec]=fixed_point(matrix);
%% Statistical setting
stat.nwalks=2;
stat.max_step=1000;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.varcut=10^(-2);
stat.vardiff=10^(-2);
dist=1;

%% Definition of the transition probability
[P, cdf]=prob_forward(H, dist);

%% Computation of the approximated solution and relative error
start=cputime;
[u_approx, var, NWALKS,reject]=MC_forward_adapt2(H, rhs, P, cdf, stat);
finish=cputime;

rel_error=norm(u-u_approx)/norm(u);

%delete(parjob)

save(strcat('../results/Adaptive_MC_forward/MC_forward_adapt2_p=', num2str(dist), '_', matrix))