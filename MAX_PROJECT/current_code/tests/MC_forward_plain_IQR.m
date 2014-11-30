addpath('../core')
addpath('../utils')

parjob=parpool('local');

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'
matrix='ifiss_convdiff';

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

n_walks=10^(4);
walkcut=10^(-6);
dist=1;
max_step=100;

[P, cdf]=prob_forward(H, dist);
%%
start=cputime;
[u_approx, var, IQR]=MC_forward_IQR(H, rhs, P, cdf, n_walks, max_step);
%[u_approx, var, X]=MC_forward(H, rhs, P, cdf, n_walks, max_step);
finish=cputime;

delete(parjob)

rel_err=norm(u-u_approx)/norm(u);

save('IQR')