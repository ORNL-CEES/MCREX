addpath('../core')
addpath('../utils')

dimen=50;
A=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
rhs=[1:50]';
u=A\rhs;

Prec=diag(diag(A));

H=eye(size(A))-Prec\A;
rhs=Prec\rhs;

dist='MAO';

[Pb, cdfb, P, cdf]=prob_adjoint(H, rhs, dist);

eps=10^(-3);
n_walks=1000;
max_step=20;
rich_it=10;

fp.u=u;
fp.H=H;
fp.rhs=rhs;
fp.precond='diag';
%% Monte Carlo Adjoint Method resolution

[sol, rel_error, var, NWALKS, iterations]=MCSA_adjoint(fp, dist, P, cdf, rich_it, n_walks, max_step, eps);

conf=0.05;

plot(sol, '*');
hold on
plot(sol-var*norminv(1-conf/2, 0, 1), 'g*');
plot(sol+var*norminv(1-conf/2, 0, 1), 'g*');
plot(u,'r*');

save(strcat('../results/MCSA_adjoint/MCSA_adjoint_test_', dist))