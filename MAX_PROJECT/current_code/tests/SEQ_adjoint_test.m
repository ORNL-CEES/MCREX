addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1', 'ifiss_convdiff'; 'shifted_laplacian_1d'; 'thermal_eq_diff'
matrix='fs_680_1';

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


Prec=diag(diag(A));

H=eye(size(A))-Prec\A;
rhs=Prec\rhs;

%% Numerical setting

numer.eps=10^(-3);
numer.rich_it=300;

%% Statistical setting

stat.nwalks=dimen;
stat.max_step=20;
stat.adapt=1;
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

[sol, rel_residual, var, VAR, DX, NWALKS, iterations]=SEQ_adjoint(fp, dist, P, cdf, numer, stat);

conf=0.05;

plot(sol, '*');
hold on
plot(sol-var*norminv(1-conf/2,0,1), 'g*');
plot(sol+var*norminv(1-conf/2,0,1), 'g*');
plot(u, 'r*');


save(strcat('../results/SEQ_adjoint/SEQ_adjoint_test_', matrix, '_p=', num2str(dist)));