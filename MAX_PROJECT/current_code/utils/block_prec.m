addpath('../core')

format long
load('SPN_mark.mat')
A=SPN_mark;

%% Block-diag preconditioning

size_block=63;
n_blocks=size(A,1)/size_block;

% D=zeros(size(A));
% 
% for i=0:n_blocks-1
%     D(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block)=A(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block);
% end

D_inv=zeros(size(A));

for i=0:n_blocks-1
    D_inv(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block)=inv(A(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block));
end

D_inv=sparse(D_inv);

B=D_inv*A;
B=sparse(B);
H=eye(size(A,1))-B;
H=sparse(H);

max(abs(eigs(H)))
P=prob_adjoint3(H,1);
MC_converge(H',P)

%% ILU preconditioning

setup.type='nofill';
setup.milu='off';

[L,U]=ilu(A,setup);

L_inv=sparse(inv(L));
U_inv=sparse(inv(U));

B=sparse(U_inv*L_inv*A);
H=eye(size(A,1))-B;
H=sparse(H);

max(abs(eigs(H)))
P=prob_adjoint3(H,0);
MC_converge(H,P)


setup.type='ilutp';
setup.milu='off';
setup.droptol = 0.2;
[L,U] = ilu(A,setup);

L_inv=sparse(inv(L));
U_inv=sparse(inv(U));

B=sparse(U_inv*L_inv*A);
H=eye(size(A,1))-B;
H=sparse(H);

max(abs(eigs(H)))
P=prob_adjoint3(H,0);
MC_converge(H,P)

%% Block Tridiagonal matrix

size_block=357;
n_blocks=size(A,1)/size_block;

D=zeros(size(A));

for i=0:n_blocks-1
    D(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block)=A(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block);
end

L=zeros(size(A));

for i=0:n_blocks-2
    L((i+1)*size_block+1:(i+2)*size_block,i*size_block+1:(i+1)*size_block)=A((i+1)*size_block+1:(i+2)*size_block,i*size_block+1:(i+1)*size_block);
end

U=zeros(size(A));

for i=0:n_blocks-2
    U(i*size_block+1:(i+1)*size_block,(i+1)*size_block+1:(i+2)*size_block)=A(i*size_block+1:(i+1)*size_block,(i+1)*size_block+1:(i+2)*size_block);
end


Delta=zeros(size(A));
Delta(1:size_block,1:size_block)=D(1:size_block,1:size_block);
aux=D(1:size_block,1:size_block);

for i=1:n_blocks-1
    Di=D(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block);
    Ei=U((i-1)*size_block+1:i*size_block,i*size_block+1:(i+1)*size_block);
    Fi=L(i*size_block+1:(i+1)*size_block,(i-1)*size_block+1:i*size_block);
    
    setup.type='nofill';
    setup.milu='off';
    %setup.droptol = 0.1;
    
    [Ltil,Util] = ilu(sparse(aux),setup);

    L_inv=sparse(inv(Ltil));
    U_inv=sparse(inv(Util));
    
    Delta(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block)=Di-Fi*U_inv*L_inv*Ei;
    
    aux= Delta(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block);
end

M=sparse((L+Delta)*Delta\(Delta+U));
B=M\A;
