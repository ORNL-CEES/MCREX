addpath('../core')
addpath('../utils')

dimen=50;
A=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
rhs=[1:50]';
u=A\rhs;

Prec=diag(diag(A));

H=eye(size(A))-Prec\A;
rhs=Prec\rhs;

eps=10^(-3);
n_walks=1000;
max_step=20;
rich_it=10;
dist='MAO';

[P, cdf]=prob_forward(H, dist);

fp.u=u;
fp.H=H;
fp.rhs=rhs;
fp.precond='diag';
%% SEQ forward method resolution

[sol, rel_err, var, NWALKS, iterations]=SEQ_forward(fp, P, cdf, rich_it, n_walks, max_step, eps);

conf=0.05;
var=var*norminv(1-conf/2, 0, 1);

count=0;
for i=1:size(u,1)
    if u(i)>(sol(i)-var(i)) && u(i)<(sol(i)+var(i))
        count=count+1;
    end
end

plot(sol, '*');
hold on
plot(sol+var, 'g*');
plot(sol-var, 'g*');
plot(u, 'r*');

save(strcat('../results/SEQ_forward/SEQ_forward_test', dist));
