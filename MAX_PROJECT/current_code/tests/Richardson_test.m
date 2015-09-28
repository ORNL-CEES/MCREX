addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'; 'sp1_ainv': 'SPN_ainv';
% 'laplacian_2d_ainv'; 'parabolic_ifiss'; 
% 'parabolic_freefemS'; 'parabolic_freefemL'; 

matrix='parabolic_ifiss'; 
[H,rhs, u, precond, Prec]=fixed_point(matrix);

%% Numerical setting

numer.eps=10^(-7);
numer.rich_it=100000;
%% Preconditioning setting

fp.u=u; %reference solution
fp.H=H; %iteration matrix
fp.rhs=rhs;
fp.precond=precond;

%% Monte Carlo Adjoint Method resolution

start=cputime;
[sol, iterations, relres]=Richardson(fp, numer);
finish=cputime;

