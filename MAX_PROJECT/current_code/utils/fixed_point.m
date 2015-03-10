function [H,rhs, precond, Prec]=fixed_point(matrix)

    addpath(strcat('../utils/model_problems/', matrix))

    if strcmp(matrix, 'simple')
        dimen=50;
        A=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
        rhs=ones(dimen,1);
        u=A\rhs;
        A=sparse(A);

        Prec=diag(diag(A));
        Prec=sparse(Prec); 
        rhs=Prec\rhs;
        H=sparse(speye(size(A))-Prec\A);
        precond='diag';

    elseif strcmp(matrix, 'SPN')
        load('A.mat');
        u=ones(size(A,1),1);
        rhs=A*u;
        size_block=69;
        n_blocks=size(A,1)/size_block;

        D_inv=zeros(size(A));

        for i=0:n_blocks-1
            D_inv(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block)=inv(A(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block));
        end

        D_inv=sparse(D_inv);

        B=D_inv*A;
        B=sparse(B);
        H=eye(size(A,1))-B;
        H=sparse(H); 
        rhs=D_inv*rhs;
        precond='block_diag';

    elseif strcmp(matrix, 'sp1')
        load('A_sp1.mat');
        A=A_sp1;
        clear A_sp1;    
        u=ones(size(A,1),1);
        rhs=A*u;
        size_block=63;
        n_blocks=size(A,1)/size_block;

        D_inv=zeros(size(A));

        for i=0:n_blocks-1
            D_inv(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block)=inv(A(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block));
        end

        D_inv=sparse(D_inv);

        B=D_inv*A;
        B=sparse(B);
        H=eye(size(A,1))-B;
        H=sparse(H); 
        rhs=D_inv*rhs;
        precond='block_diag';

    elseif strcmp(matrix, 'sp1_shift') 
         load('A.mat')
         u=ones(size(A,1),1);
         rhs=A*u;
         Prec=sparse(diag(diag(A)));
         H=speye((size(A,1)),(size(A,1)))-A*sparse(inv(Prec));
         precond='diag';

    elseif strcmp(matrix, 'sp3_shift') 
         load('A.mat')
         u=ones(size(A,1),1);
         rhs=A*u;
         Prec=sparse(diag(diag(A)));
         H=speye((size(A,1)),(size(A,1)))-A*sparse(inv(Prec));
         precond='diag';


    elseif strcmp(matrix, 'SPN_shift')
         load('A.mat')
         u=ones(size(A,1),1);
         rhs=A*u;
         Prec=sparse(diag(diag(A)));
         H=speye((size(A,1)),(size(A,1)))-A*sparse(inv(Prec));
         precond='diag';

    elseif strcmp(matrix, 'sp5_shift')
         load('A.mat')
         u=ones(size(A,1),1);
         rhs=A*u;
         Prec=sparse(diag(diag(A)));
         H=speye((size(A,1)),(size(A,1)))-A*sparse(inv(Prec));
         precond='diag';
    elseif strcmp(matrix, 'sp1_ainv') || strcmp(matrix, 'SPN_ainv')
        load('A.mat')
        load('prec_z')
        Z_t=spconvert(prec_z);
        Z_t=sparse(Z_t);
        load('prec_w')
        W=spconvert(prec_w);
        W=sparse(W);
        nn=size(Z_t,1);
        W=W(1:nn,1:nn);
        Z_t=Z_t(1:nn,1:nn);

        D_inv=sparse(diag(diag(Z_t)));
        Z=Z_t'*inv(D_inv);
        C=sparse(Z * D_inv * W);
        H = sparse(speye(nn) - A*C);
        
        u=ones(size(A,1),1);
        rhs=A*u; 
        Prec=C;
        precond='ainv';

    else
         [A, dimen, ~, ~] = mmread('A.mtx');
         rhs=mmread('b.mtx');
         u=mmread('x.mtx');
         A=sparse(A);
         Prec=diag(diag(A));
         Prec=sparse(Prec);
         H=sparse(speye(size(A))-Prec\A);
         rhs=Prec\rhs;
         precond='diag';
    end


return