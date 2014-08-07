addpath('../core')
addpath('../utils')

dimen=500;
A=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
rhs=[1:500]';
u=A\rhs;

Prec=diag(diag(A));

H=eye(size(A))-Prec\A;
rhs=Prec\rhs;

max_step=100;
walkcut=10^(-6);
dist_un=0;
dist_mao=1;

[P_un, cdf_un]=prob_forward(H, dist_un);
[P_mao, cdf_mao]=prob_forward(H, dist_mao);

%%
n_walks=[10^1 10^2 10^3 10^4 10^5];

[sol_un, var_un, err_un, time_un]=MC_forward_error2(u, H, rhs, P_un, cdf_un, n_walks, max_step, walkcut);

rel_error_un=[];
for i=1:length(n_walks)
    rel_error_un=[rel_error_un sqrt(sum((u-sol_un(:,i)).^2))/sqrt(sum((u.^2)))];
end

[sol_mao, var_mao, err_mao, time_mao]=MC_forward_error2(u, H, rhs, P_mao, cdf_mao, n_walks, max_step, walkcut);

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

FoM_un=1./((norm_var_un.^2).*time_un);
FoM_mao=1./((norm_var_mao.^2).*time_mao);

hold off
loglog(n_walks, sqrt(1./n_walks), 'r');
hold on
loglog(n_walks, rel_error_un, '-ok');
hold on
loglog(n_walks, rel_error_mao, '-ob');

hold off
loglog(n_walks, sqrt(1./n_walks), 'r');
hold on
loglog(n_walks, norm_var_un, '-ok');
hold on
loglog(n_walks, norm_var_mao, '-ob');

hold off
semilogx(n_walks, FoM_un, '-or');
hold on
semilogx(n_walks, FoM_mao, '-ob');

save(strcat('../results/MC_forward_compare/MC_forward_plain'))
%save(strcat('MC_forward_plain_', dist))