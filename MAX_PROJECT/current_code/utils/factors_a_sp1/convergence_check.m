addpath('../../core');

load('A_sp1.mat')
load('z')
Z_t=spconvert(z);
Z_t=sparse(Z_t);
load('w')
W=spconvert(w);
W=sparse(W);
D_inv=sparse(diag(diag(Z_t)));
Z_t=sparse( Z_t - D_inv + sparse(diag(ones(18207,1))) );
Z=Z_t';
C=sparse(Z * D_inv * W);

%% Study of positive definitiness
B=A_sp1+A_sp1';
p=symamd(B);
Chol=chol(B(p,p));


%% Study of the spectral radius
H = sparse(speye(18207) - C*A_sp1);

b=zeros(size(H,1),1);
[P] = prob_adjoint3(H, 1);
[answer]=MC_converge(H,P)
