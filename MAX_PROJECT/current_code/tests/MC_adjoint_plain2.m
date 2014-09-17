addpath('../core')
addpath('../utils')

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'
matrix='jpwh_991';

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
dist2=1;

[Pb, cdfb, P, cdf]=prob_adjoint2(H, rhs, dist1, dist2);
%%
n_walks=[10^1 10^2 10^3 10^4];

max_step=1000;

walkcut1=1;
[sol1, var1, tally1, time1]=MC_adjoint_error2(H, rhs, P, cdf, Pb, cdfb, n_walks, max_step, walkcut1);

walkcut2=10^(-3);
[sol2, var2, tally2, time2]=MC_adjoint_error2(H, rhs, P, cdf, Pb, cdfb, n_walks, max_step, walkcut2);

walkcut3=10^(-6);
[sol3, var3, tally3, time3]=MC_adjoint_error2(H, rhs, P, cdf, Pb, cdfb, n_walks, max_step, walkcut3);

walkcut4=10^(-9);
[sol4, var4, tally4, time4]=MC_adjoint_error2(H, rhs, P, cdf, Pb, cdfb, n_walks, max_step, walkcut4);
    
rel_error1=[];
for i=1:length(n_walks)
    rel_error1=[rel_error1 sqrt(sum((u-sol1(:,i)).^2))/sqrt(sum((u.^2)))];
end

rel_error2=[];
for i=1:length(n_walks)
    rel_error2=[rel_error2 sqrt(sum((u-sol2(:,i)).^2))/sqrt(sum((u.^2)))];
end

rel_error3=[];
for i=1:length(n_walks)
    rel_error3=[rel_error3 sqrt(sum((u-sol3(:,i)).^2))/sqrt(sum((u.^2)))];
end

rel_error4=[];
for i=1:length(n_walks)
    rel_error4=[rel_error4 sqrt(sum((u-sol4(:,i)).^2))/sqrt(sum((u.^2)))];
end


loglog(n_walks, rel_error1, '-or');
hold on
loglog(n_walks, 1./sqrt(n_walks), 'k');
loglog(n_walks, rel_error2, '-og');
loglog(n_walks, rel_error3,'-ob');
loglog(n_walks, rel_error4,'-oc')
legend('W_c=1', '1/sqrt(N)', 'W_c=0.1', 'W_c=0.01', 'W_c=0.001')
title('RELATIVE ERROR - SHIFTED LAPLACIAN')
xlabel('Nb. Random Walks');
ylabel('Rel. Error')


save(strcat('../results/MC_adjoint_plain2/MC_adjoint_plain2_p1=', num2str(dist1), '_p2=', num2str(dist2), '_', matrix));
