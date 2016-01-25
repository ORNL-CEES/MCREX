addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'; 'sp1_ainv': 'SPN_ainv';
% 'laplacian_2d_ainv'; 'parabolic_ifiss'; 
% 'parabolic_freefemS'; 'parabolic_freefemL'; 
% 'parabolic_freefemS_diag'; 'parabolic_freefemL_diag'; 

matrix='laplacian_2d';
[H,rhs, u, precond, Prec]=fixed_point(matrix);

%% Statistical setting
stat.nwalks=10;
stat.max_step=1000;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.varcut=10^(-3);
stat.vardiff=10^(-2);
dist=1;

%% Definition of the transition probability
[Pb, cdfb, P, cdf]=prob_adjoint(H, rhs, dist);

%% Computation of the approximated solution and relative error
start=cputime;
[u_approx, var, tot_nwalks, tally, res]=MC_adjoint_adapt(H, rhs, P, cdf, Pb, cdfb, stat);
finish=cputime;

rel_error=norm(u-u_approx)/norm(u);

if strcmp(matrix, 'sp1_shift') || strcmp(matrix, 'sp3_shift') || strcmp(matrix, 'SPN_shift') || strcmp(matrix, 'sp5_shift') || ...
     strcmp(matrix, 'parabolic_freefemL_diag') || strcmp(matrix, 'parabolic_freefemS_diag')    
    u_approx=Prec\u_approx;
    
elseif strcmp(matrix, 'sp1_ainv') || strcmp(matrix, 'SPN_ainv')  || strcmp(matrix, 'parabolic_ifiss')
    u_approx=Prec*u_approx;
end

%delete(parjob)

save(strcat('../results/Adaptive_MC_adjoint/MC_adjoint_adapt2_p=', num2str(dist), '_', matrix))
