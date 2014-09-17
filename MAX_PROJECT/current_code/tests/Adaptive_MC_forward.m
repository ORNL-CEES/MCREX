addpath('../core')
addpath('../utils')

parjob=parpool('local');

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'
matrix='laplacian_2d';

if ~strcmp(matrix, 'simple')
    addpath(strcat('../utils/model_problems/', matrix));
end


if strcmp(matrix, 'simple')
    dimen=500;
    A=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
    rhs=[1:dimen]';
    u=A\rhs;

else
     [A, dimen, ~, ~] = mmread('A.mtx');
     rhs=mmread('b.mtx');
     u=mmread('x.mtx');
end

Prec=diag(diag(A));

H=eye(size(A))-Prec\A;
rhs=Prec\rhs;

%% Statistical setting
stat.nwalks=2;
stat.max_step=1000;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.varcut=10^(-1);
stat.vardiff=10^(-1);
dist=1;

%% Definition of the transition probability
[P, cdf]=prob_forward(H, dist);

%% Computation of the approximated solution and relative error
start=cputime;
[u_approx, var, NWALKS, reject]=MC_forward_adapt2(H, rhs, P, cdf, stat);
finish=cputime;

rel_error=sqrt(sum((u-u_approx).^2))/sqrt(sum((u.^2)));

delete(parjob)

save(strcat('../results/Adaptive_MC_forward/MC_forward_adapt2_p=', num2str(dist), '_', matrix))