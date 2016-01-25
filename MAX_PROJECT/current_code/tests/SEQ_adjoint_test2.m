addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'; 'sp1_ainv': 'SPN_ainv';

matrix='thermal_eq_diff';

[H,rhs, u, precond, Prec]=fixed_point(matrix);

%% Numerical setting

numer.eps=10^(-3);
numer.rich_it=300;

%% Statistical setting
stat.nwalks=2;
stat.max_step=20;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.nchecks=2;
stat.varcut=0.5;
stat.vardiff=10^(-1);
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
[sol, rel_residual, ~, ~, ~, NWALKS, tally, iterations, reject]=SEQ_adjoint(fp, dist, P, cdf, numer, stat);
finish=cputime;


% for i=1:iterations
%     figure()
%     norm_var=[];
%     for j=1:size(VAR{i},2)
%         norm_var=[norm_var norm(VAR{i}(:,j))];
%     end
%     loglog(norm_var, '-o')
% end
% 
% for i=1:iterations
%     figure()
%     loglog(RES{i}, '-o')
%     %axis([1 length(norm_res) 10^(-1) 10^0]);
% end

if strcmp(matrix, 'sp1_shift') || strcmp(matrix, 'sp3_shift') || strcmp(matrix, 'SPN_shift') || strcmp(matrix, 'sp5_shift') || ...
     strcmp(matrix, 'parabolic_freefemL_diag') || strcmp(matrix, 'parabolic_freefemS_diag')    
    sol=Prec\sol;
    
elseif strcmp(matrix, 'sp1_ainv') || strcmp(matrix, 'SPN_ainv')  || strcmp(matrix, 'parabolic_ifiss')
    sol=Prec*sol;
end

save(strcat('../results/MCSA_adjoint2/MCSA_adjoint_test2_', matrix, '_p=', num2str(dist)))