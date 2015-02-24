addpath('../core')
addpath('../utils')


% 'jpwh_991'; 'fs_680_1'; 'ifiss_convdiff'; 'shifted_laplacian_1d';
% 'thermal_eq_diff'; 'laplacian_2d'; 'SPN'; 'sp1'; 'SPN_shift';
% 'sp1_shift'; 'sp5_shift'; 'sp3_shift'
matrix= 'ifiss_convdiff';
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
    
elseif strcmp(matrix, 'sp1_shift') 
     load('A.mat')
     u=ones(size(A,1),1);
     rhs=A*u;
     Prec=sparse(diag(diag(A)));
     H=speye((size(A,1)),(size(A,1)))-Prec\A;
     rhs=Prec\rhs;
    
elseif strcmp(matrix, 'sp3_shift') 
     load('A.mat')
     u=ones(size(A,1),1);
     rhs=A*u;
     Prec=sparse(diag(diag(A)));
     H=speye((size(A,1)),(size(A,1)))-Prec\A;
     rhs=Prec\rhs;
          
     
elseif strcmp(matrix, 'SPN_shift')
     load('A.mat')
     u=ones(size(A,1),1);
     rhs=A*u;
     Prec=sparse(diag(diag(A)));
     H=speye((size(A,1)),(size(A,1)))-Prec\A;
     rhs=Prec\rhs;  
     
elseif strcmp(matrix, 'sp5_shift')
     load('A.mat')
     u=ones(size(A,1),1);
     rhs=A*u;
     Prec=sparse(diag(diag(A)));
     H=speye((size(A,1)),(size(A,1)))-Prec\A;
     rhs=Prec\rhs;        
    
else
     [A, dimen, ~, ~] = mmread('A.mtx');
     rhs=mmread('b.mtx');
     u=mmread('x.mtx');
     A=sparse(A);
     Prec=diag(diag(A));
     Prec=sparse(Prec);
     H=sparse(speye(size(A))-Prec\A);
     rhs=Prec\rhs;
end

walkcut=10^(-6);
dist=1;

[P, cdf]=prob_forward(H, dist);

max_step=1000;

%%
n_walks=[10^1 10^2 10^3 10^4];

start=cputime;
[u_approx, var, err]=MC_forward_error(u, H, rhs, P, cdf, n_walks, max_step);
finish=cputime;

delete(parjob)

rel_error=[];
for i=1:length(n_walks)
    rel_error=[rel_error sqrt(sum((u-u_approx(:,i)).^2))/sqrt(sum((u.^2)))];
end



loglog(n_walks, sqrt(1./n_walks), 'r')
hold on
for i=1:size(rhs,1)
    loglog(n_walks, err(:,i), 'b')
    hold on
end
loglog(n_walks, rel_error, '-*g', 'Linewidth', 3)

hold off
figure()
loglog(n_walks, sqrt(1./n_walks), 'r');
hold on
loglog(n_walks, rel_error, '-o');

conf=0.05;

count=0;

for i=1:size(u,1)
    if u(i)>u_approx(i)-var(i)*norminv(1-conf/2, 0, 1) && u(i)<u_approx(i)+var(i)*norminv(1-conf/2, 0, 1)
        count=count+1;
    end
end

amplitude=0;
for i=1:size(u,1)
    amplitude=amplitude+2*var(i)*norminv(1-conf/2, 0, 1);
end
amplitude=amplitude/norm(u,2);


hold off
plot(u_approx(:,end), 'k*')
hold on
plot(u_approx(:,end)-var(:,end)*norminv(1-conf/2, 0,1), 'g*')
plot(u_approx(:,end)+var(:,end)*norminv(1-conf/2, 0,1), 'g*')
plot(u, 'r*')


save(strcat('../results/MC_forward_plain/MC_forward_plain_p=', num2str(dist), '_', matrix))
%save(strcat('MC_forward_plain_', dist))