load('A_sp1_pre2.mat')
load('A_sp1.mat')

[Z,W,D_inv]=split_michele(A_sp1_pre2);
Prec=Z*D_inv*W;

atv=@(x) A_sp1*x;

x0=zeros(size(A_sp1,1),1);
b=A_sp1*ones(size(A_sp1,1),1);
tol=10^(-9);
max_iter=200;
params=[tol, max_iter];
[x, error, total_iters] = gmres2(x0, b, atv, params, Prec);
