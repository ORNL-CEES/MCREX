addpath('../core')
addpath('../utils')

parjob=parpool('local');

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'
matrix='fs_680_1';

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

walkcut=10^(-6);
dist=1;

[P, cdf]=prob_forward(H, dist);

max_step=10;

%%
n_walks=[10 100 1000];

start=cputime;
[u_approx, var, err]=MC_forward_error(u, H, rhs, P, cdf, n_walks, max_step);
finish=cputime;

delete(parjob)

rel_error=[];
for i=1:length(n_walks)
    rel_error=[rel_error sqrt(sum((u-u_approx(:,i)).^2))/sqrt(sum((u.^2)))];
end



loglog(n_walks, sqrt(1./n_walks), 'r')
hold on
for i=1:size(rhs,1)
    loglog(n_walks, err(:,i), 'b')
    hold on
end
loglog(n_walks, rel_error, '-*g', 'Linewidth', 3)

hold off
figure()
loglog(n_walks, sqrt(1./n_walks), 'r');
hold on
loglog(n_walks, rel_error, '-o');

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
plot(u_approx(:,end), 'k*')
hold on
plot(u_approx(:,end)-var(:,end)*norminv(1-conf/2, 0,1), 'g*')
plot(u_approx(:,end)+var(:,end)*norminv(1-conf/2, 0,1), 'g*')
plot(u, 'r*')


save(strcat('../results/MC_forward_plain/MC_forward_plain_p=', num2str(dist), '_', matrix))
%save(strcat('MC_forward_plain_', dist))