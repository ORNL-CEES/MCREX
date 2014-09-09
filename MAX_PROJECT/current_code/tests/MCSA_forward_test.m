addpath('../core')
addpath('../utils')

parjob=parpool('local');

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff';  'laplacian_2d'
matrix='laplacian_2d';

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

%% NUmerical setting

numer.eps=10^(-3);
numer.rich_it=300;

%% Statistical setting

stat.nwalks=10;
stat.max_step=1000;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-18);
stat.nchecks=1;
stat.varcut=0.1;
dist=1;

%% Definition of the transitional probability
[P, cdf]=prob_forward(H, dist);

%% Preconditioning setting

fp.u=u; %reference solution
fp.H=H; %iteration matrix
fp.rhs=rhs;
fp.precond='diag';
%% MCSA forward method resolution

start=cputime;
[sol, rel_residual, var, VAR, DX, NWALKS, tally, iterations, reject]=MCSA_forward(fp, P, cdf, numer, stat);
finish=cputime;

delete(parjob);

conf=0.05;

count=0;

for i=1:size(u,1)
    if u(i)>sol(i)-var(i)*norminv(1-conf/2, 0, 1) && u(i)<sol(i)+var(i)*norminv(1-conf/2, 0, 1)
        count=count+1;
    end
end

plot(sol, '*');
hold on
plot(sol-VAR(:,end)*norminv(1-conf/2, 0, 1), 'g*');
plot(sol+VAR(:,end)*norminv(1-conf/2, 0, 1), 'g*');
plot(u, 'r*')

hold off
bar(NWALKS)

save(strcat('../results/MCSA_forward/MCSA_forward_test_', matrix, '_p=', num2str(dist)))