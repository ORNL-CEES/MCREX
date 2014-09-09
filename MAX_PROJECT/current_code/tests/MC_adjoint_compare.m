addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'
matrix='simple';

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

eps=10^(-3);
dist1=1;
dist_un=0;
dist_mao=1;
walkcut=10^(-6);
max_step=1000;

n_walks=[10^1 10^2 10^3 10^4]; %10^5 10^6];

[Pb, cdfb, P_un, cdf_un]=prob_adjoint2(H, rhs, dist1, dist_un);
[~, ~, P_mao, cdf_mao]=prob_adjoint2(H, rhs, dist1, dist_mao);
%%
[sol_un, var_un, tally_un, time_un]=MC_adjoint_error(H, rhs, P_un, cdf_un, Pb, cdfb, n_walks, max_step);

rel_error_un=[];
for i=1:length(n_walks)
    rel_error_un=[rel_error_un sqrt(sum((u-sol_un(:,i)).^2))/sqrt(sum((u.^2)))];
end

%%
[sol_mao, var_mao, tally_mao, time_mao]=MC_adjoint_error(H, rhs, P_mao, cdf_mao, Pb, cdfb, n_walks, max_step);

rel_error_mao=[];
for i=1:length(n_walks)
    rel_error_mao=[rel_error_mao sqrt(sum((u-sol_mao(:,i)).^2))/sqrt(sum((u.^2)))];
end

norm_var_un=[];
for i=1:length(n_walks)
    norm_var_un=[norm_var_un norm(var_un(:,i).^2,1)];
end

norm_var_mao=[];
for i=1:length(n_walks)
    norm_var_mao=[norm_var_mao norm(var_mao(:,i).^2,1)];
end
%% 
FoM_un=1./((norm_var_un).*time_un);
FoM_mao=1./((norm_var_mao).*time_mao);

figure()   
loglog(n_walks, rel_error_un, '-or');
hold on
loglog(n_walks, 1./sqrt(n_walks), 'k')
hold on
loglog(n_walks, rel_error_mao, '-ob');
title('ADJOINT MONTE CARLO - ERROR BEHAVIOUR');
xlabel('Nb. random walks');
ylabel('Rel. Error');

hold off
figure()
loglog(n_walks, sqrt(1./n_walks), 'r');
hold on
loglog(n_walks, norm_var_un, '-ok');
hold on
loglog(n_walks, norm_var_mao, '-ob');

hold off
figure()
semilogx(n_walks, FoM_un, '-or');
hold on
semilogx(n_walks, FoM_mao, '-ob');
title('FIGURE OF MERIT')
xlabel('Nb. random walks');
ylabel('FoM');

save(strcat('../results/MC_adjoint_compare/MC_adjoint_compare_', matrix));
%   save(strcat('MC_adjoint_plain_', dist)); 
