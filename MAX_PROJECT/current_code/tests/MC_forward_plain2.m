addpath('../core')
addpath('../utils')

parjob=parpool('local');

% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'
matrix='laplacian_2d';

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

dist=1;

[P, cdf]=prob_forward(H, dist);

max_step=1000;

%%
n_walks=[10 100 1000 10000];

walkcut1=10^(0);
start1=cputime;
[sol1, var1, err1]=MC_forward_error2(u, H, rhs, P, cdf, n_walks, max_step, walkcut1);
finish1=cputime;

rel_error1=[];
for i=1:length(n_walks)
    rel_error1=[rel_error1 sqrt(sum((u-sol1(:,i)).^2))/sqrt(sum((u.^2)))];
end


walkcut2=10^(-3);
start2=cputime;
[sol2, var2, err2]=MC_forward_error2(u, H, rhs, P, cdf, n_walks, max_step, walkcut2);
finish2=cputime;

rel_error2=[];
for i=1:length(n_walks)
    rel_error2=[rel_error2 sqrt(sum((u-sol2(:,i)).^2))/sqrt(sum((u.^2)))];
end

walkcut3=10^(-6);
start3=cputime;
[sol3, var3, err3]=MC_forward_error2(u, H, rhs, P, cdf, n_walks, max_step, walkcut3);
finish3=cputime;

rel_error3=[];
for i=1:length(n_walks)
    rel_error3=[rel_error3 sqrt(sum((u-sol3(:,i)).^2))/sqrt(sum((u.^2)))];
end

walkcut4=10^(-9);
start4=cputime;
[sol4, var4, err4]=MC_forward_error2(u, H, rhs, P, cdf, n_walks, max_step, walkcut4);
finish4=cputime;

rel_error4=[];
for i=1:length(n_walks)
    rel_error4=[rel_error4 sqrt(sum((u-sol4(:,i)).^2))/sqrt(sum((u.^2)))];
end

delete(parjob)

hold off
loglog(n_walks, sqrt(1./n_walks), 'k');
hold on
loglog(n_walks, rel_error1, '-or');
loglog(n_walks, rel_error2, '-ob');
loglog(n_walks, rel_error3, '-og');
loglog(n_walks, rel_error4, '-oc');
title('RELATIVE ERROR - FORWARD MONTE CARLO - SHIFTED LAPLACIAN');
xlabel('Nb. Random Walks');
ylabel('Relative Error');
legend('1/sqrt(N)', 'W_c=1', 'W_c=0.01', 'W_c=0.001', 'W_c=10^(-6)');

save(strcat('../results/MC_forward_plain2/MC_forward_plain2_p=', num2str(dist), '_', matrix))
