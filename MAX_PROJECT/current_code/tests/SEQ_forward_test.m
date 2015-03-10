addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'; 'sp1_ainv': 'SPN_ainv';

matrix='thermal_eq_diff';

[H,rhs, precond, Prec]=fixed_point(matrix);

%% Numerical setting

numer.eps=10^(-3);
numer.rich_it=300;

%% Statistical setting

stat.nwalks=2;
stat.max_step=20;
stat.adapt_walks=1;
stat.adapt_cutoff=1;
stat.walkcut=10^(-6);
stat.varcut=0.5;
dist=1;

%% Definition of the transitional probability
[P, cdf]=prob_forward(H, dist);

%% Preconditioning setting

fp.u=u; %reference solution
fp.H=H; %iteration matrix
fp.rhs=rhs;
fp.precond='diag';
%% SEQ forward method resolution

[sol, rel_residual, var, VAR, DX, NWALKS, tally, iterations]=SEQ_forward(fp, P, cdf, numer, stat);

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

save(strcat('../results/SEQ_forward/SEQ_forward_test_', matrix, '_p=', num2str(dist)));
