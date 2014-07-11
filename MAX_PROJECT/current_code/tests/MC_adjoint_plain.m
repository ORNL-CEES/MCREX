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
dist1=1;
dist2=1;

[Pb, cdfb, P, cdf]=prob_adjoint2(H, rhs, dist1, dist2);
%%
rel_errs=[];
NWALKS=[];
VAR=[];
SOL=[];

for i=1:2
    X=[];
    n_walks=10^i;
    NWALKS=[NWALKS n_walks];
    max_step=100;
    [u_approx, var, X, tally]=MC_adjoint(X,H,rhs,P,cdf,Pb,cdfb, n_walks, max_step);
    rel_error=sqrt(sum((u-u_approx).^2))/sqrt(sum((u).^2));
    rel_errs=[rel_errs, rel_error];
    SOL=[SOL u_approx];
    VAR=[VAR var];
end
    
loglog(NWALKS, rel_errs, '-or');
hold on
loglog(NWALKS, 1./sqrt(NWALKS), 'k')

conf=0.05;

count=0;

for i=1:size(u,1)
    if u(i)>u_approx(i)-var(i)*norminv(1-conf/2, 0, 1) && u(i)<u_approx(i)+var(i)*norminv(1-conf/2, 0, 1)
        count=count+1;
    end
end

amplitude=0;
for i=1:size(u,1)
    amplitude=amplitude+2*var(i)*norminv(1-conf/2, 0, 1);
end
amplitude=amplitude/norm(u,2);

hold off
plot(u_approx, 'k*')
hold on
plot(u_approx-var*norminv(1-conf/2, 0,1), 'g*')
plot(u_approx+var*norminv(1-conf/2, 0,1), 'g*')
plot(u, 'r*')


save(strcat('../results/MC_adjoint_plain/MC_adjoint_plain_p1=', num2str(dist1), '_p2=', num2str(dist2)));
%   save(strcat('MC_adjoint_plain_', dist)); 
