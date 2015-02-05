format long
load('A_sp1.mat')
load('z')
Z_t=spconvert(z);
Z_t=sparse(Z_t);
load('w')
W=spconvert(w);
W=sparse(W);
nn=18207;
W=W(1:nn,1:nn);
Z_t=Z_t(1:nn,1:nn);
A_sp1=A_sp1(1:nn,1:nn);

D_inv=sparse(diag(diag(Z_t)));
Z=Z_t'*inv(D_inv);
C=sparse(Z * D_inv * W);

%% Study of positive definitiness
B=A_sp1+A_sp1';
p=symamd(B);
Chol=chol(B(p,p));


%% Study of the spectral radius
H = sparse(speye(nn) - C*A_sp1);

display('\n Spectral radius of original iteration matrix --- left preconditioning:\n');
max(abs(eigs(H)))

%% tests about iterative methods
rhs=A_sp1*ones(nn,1);
% check with the bicgstab
max_itbcg=100;
max_itgm=100;
tol=10^(-8);

% defiition of the left preconditioner
Prec_A=C*A_sp1;
Prec_rhs=C*rhs;

[x_bicg,flag_bicg,relres_bicg,iter_bicg,resvec_bicg] = bicgstab(Prec_A,Prec_rhs,10^(-8)*norm(Prec_rhs), max_itbcg);
display('Relative error for BiCGASTAB --- left preconditioning:')
norm((ones(nn,1))-x_bicg)/norm(ones(nn,1))
semilogy(resvec_bicg/norm(Prec_rhs))
title('BiCGASTAB --- left preconditioning: trend of the relative redisual for all the half iterations computed');

[x_gm,flag_gm,relres_gm,iter_gm, resvec_gm] = gmres(Prec_A,Prec_rhs,[],10^(-8)*norm(Prec_rhs), max_itgm);
display('Relative error for GMRES --- left preconditioning:')
norm((ones(nn,1))-x_gm)/norm(ones(nn,1))
figure()
semilogy(resvec_gm/norm(Prec_rhs))
title('GMRES --- left preconditioning: trend of the relative redisual for all the half iterations computed');