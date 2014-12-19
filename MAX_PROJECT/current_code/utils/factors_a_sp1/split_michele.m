format long
load('A_sp1.mat')
load('prec_z')
Z_t=spconvert(prec_z);
Z_t=sparse(Z_t);
load('prec_w')
W=spconvert(prec_w);
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
H = sparse(speye(18207) - C*A_sp1);

display('\n Spectral radius of original iteration matrix --- left preconditioning:\n');
max(abs(eigs(H)))

C_tilde=filtering(C, 5*10^(-1));
H_tilde = sparse(speye(18207) - C*A_sp1);

display('\n Spectral radius of filtered iteration matrix --- left preconditioning:\n');
max(abs(eigs(H_tilde)))

Z2=filtering(Z, 10^(-1));
W2=filtering(W, 10^(-1));
C2=sparse(Z2 * D_inv * W2);
H2 = sparse(speye(18207) - C2*A_sp1);

display('\n Spectral radius of iteration matrix after filtering each factor --- left preconditioning:\n');
max(abs(eigs(H2)))

display('Spectral radius of iteration matrix with two level filtering --- left preconditioning:\n');
C3=sparse(filtering( Z2 * D_inv * W2, 10^(-1)));
H3 = sparse(speye(18207) - C3*A_sp1);
max(abs(eigs(H3)))

H4 = speye(18207) - sparse(D_inv*W*A_sp1*Z);
display('\n Spectral radius of original iteration matrix --- double sided preconditioning:\n');
max(abs(eigs(H4)))


H5 = sparse( speye(18207) - D_inv*W2*A_sp1*Z2 );
display('\n Spectral radius of iteration matrix after filtering each factor--- double sided preconditioning:\n');
max(abs(eigs(H5)))


%% tests about iterative methods

rhs=A_sp1*ones(18207,1);
% check with the bicgstab
max_itbcg=100;
max_itgm=100;
tol=10^(-8);

% defiition of the left preconditioner
Prec_A=C*A_sp1;
Prec_rhs=C*rhs;

[x_bicg,flag_bicg,relres_bicg,iter_bicg,resvec_bicg] = bicgstab(Prec_A,Prec_rhs,10^(-8)*norm(Prec_rhs), max_itbcg);
display('Relative error for BiCGSTAB --- left preconditioning:')
norm((ones(18207,1))-x_bicg)/norm(ones(18207,1))
semilogy(resvec_bicg/norm(Prec_rhs))
title('BiCGASTAB --- left preconditioning: trend of the relative redisual for all the half iterations computed');

[x_gm,flag_gm,relres_gm,iter_gm, resvec_gm] = gmres(Prec_A,Prec_rhs,[],10^(-8)*norm(Prec_rhs), max_itgm);
display('Relative error for GMRES --- left preconditioning:')
norm((ones(18207,1))-x_gm)/norm(ones(18207,1))
figure()
semilogy(resvec_gm/norm(Prec_rhs))
title('GMRES --- left preconditioning: trend of the relative redisual for all the half iterations computed');

%definition of the right preconditioner
Prec_A2=A_sp1*C;

[y_bicg2,flag_bicg2,relres_bicg2,iter_bicg2,resvec_bicg2] = bicgstab(Prec_A2,rhs,10^(-8)*norm(rhs), max_itbcg);
display('Relative error for BiCGSTAB --- right preconditioning:')
x_bicg2=C*y_bicg2;
norm((ones(18207,1))-x_bicg2)/norm(ones(18207,1))
figure()
semilogy(resvec_bicg2/norm(rhs))
title('BiCGASTAB --- right preconditioning: trend of the relative redisual for all the half iterations computed');

[y_gm2,flag_gm2,relres_gm2,iter_gm2, resvec_gm2] = gmres(Prec_A2,rhs,[],10^(-8)*norm(rhs), max_itgm);
display('Relative error for GMRES:')
x_gm2=C*y_gm2;
norm((ones(18207,1))-x_gm2)/norm(ones(18207,1))
figure()
semilogy(resvec_gm2/norm(rhs))
title('GMRES --- right preconditioning: trend of the relative redisual for all the half iterations computed');



Prec_A3=D_inv*W*A_sp1*Z;
Prec_rhs3=D_inv*W*rhs;

[y_bicg3,flag_bicg3,relres_bicg3,iter_bicg3,resvec_bicg3] = bicgstab(Prec_A3,Prec_rhs3,10^(-8)*norm(Prec_rhs3), max_itbcg);
display('Relative error for BiCGSTAB --- double sided preconditioning:')
x_bicg3=Z*y_bicg3;
norm((ones(18207,1))-x_bicg3)/norm(ones(18207,1))
figure()
semilogy(resvec_bicg3/norm(Prec_rhs3))
title('BiCGASTAB --- double sided preconditioning: trend of the relative redisual for all the half iterations computed');

[y_gm3,flag_gm3,relres_gm3,iter_gm3, resvec_gm3] = gmres(Prec_A3,Prec_rhs3,[],10^(-8)*norm(Prec_rhs3), max_itgm);
display('Relative error for GMRES --- double sided preconditioning:')
x_gm3=Z*y_gm3;
norm((ones(18207,1))-x_gm3)/norm(ones(18207,1))
figure()
semilogy(resvec_gm3/norm(Prec_rhs3))
title('GMRES --- double sided preconditioning: trend of the relative redisual for all the half iterations computed');

%% diagonal preconditioner
Jacobi=sparse(diag(1./diag(A_sp1)));

Prec_A4=sparse(Jacobi*A_sp1);
Prec_rhs4=Jacobi*rhs;
max_itdiag=100;

[y_bicg4,flag_bicg4,relres_bicg4,iter_bicg4,resvec_bicg4] = bicgstab(Prec_A4,Prec_rhs4,10^(-8)*norm(Prec_rhs4), max_itdiag);



save('SALVAMI');