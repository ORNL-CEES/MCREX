addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'; 'sp1_ainv': 'SPN_ainv';
% 'laplacian_2d_ainv'; 'parabolic_ifiss'; 
% 'parabolic_freefemS'; 'parabolic_freefemL'; 
% 'parabolic_freefemS_diag'; 'parabolic_freefemL_diag'; 

matrix='jpwh_991';
[H,rhs, u, precond, Prec]=fixed_point(matrix);

%% Numerical setting

numer.eps=10^(-7);
numer.rich_it=100000;

%% Statistical setting

stat.nwalks=2;
stat.max_step=10;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.nchecks=1;
stat.varcut=0.01;
dist=1;

%% Definition of the transitional probability
[P, cdf]=prob_forward(H, dist);

%% Preconditioning setting

fp.u=u; %reference solution
fp.H=H; %iteration matrix
fp.rhs=rhs;
fp.precond=precond;

%% MCSA forward method resolution

start=cputime;
%parjob=parpool('local');
%[sol, rel_residual, VAR, DX, NWALKS, tally, iterations, reject]=MCSA_forward(fp, P, cdf, numer, stat);
[sol, rel_residual, ~, ~, ~, NWALKS, ~, iterations, reject]=MCSA_forward(fp, P, cdf, numer, stat);
%delete(parjob);
finish=cputime;

if strcmp(matrix, 'sp1_shift') || strcmp(matrix, 'sp3_shift') || strcmp(matrix, 'SPN_shift') || strcmp(matrix, 'sp5_shift') || ...
     strcmp(matrix, 'parabolic_freefemL_diag') || strcmp(matrix, 'parabolic_freefemS_diag')    
    sol=Prec\sol;
    
elseif strcmp(matrix, 'sp1_ainv') || strcmp(matrix, 'SPN_ainv')  || strcmp(matrix, 'parabolic_ifiss')
    sol=Prec*sol;
end

bar(NWALKS)

save(strcat('../results/MCSA_forward/MCSA_forward_test_', matrix, '_p=', num2str(dist)))