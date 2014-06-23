function [fp]=iteration_matrix(precond, D, u, rhs, G)

fp.precond=precond;
fp.D=D;
MATRICES=struct;
SET=struct;

if strcmp(precond, 'diag') 
    [~, fp]=Prec(D, precond, fp, SET);
    fp.H=eye(size(D))-(fp.prec)\D;
    fp.rhs=(fp.prec)\rhs;
    fp.u=u;

elseif strcmp(precond, 'triblock')
    %reordering
    fp.Per=redblack(G);
   
    [MATRICES, fp]=Prec(D, precond, fp, MATRICES);
    
    fp.M=MATRICES.M;
    fp.N=MATRICES.N;
       
    %building of the iteration matrix
    Pr=fp.M;
    H=eye(size(D,1))-Pr\D;
    rhs=Pr\rhs;
    Num=size(H,1);
    H(1:Num/2, 1:Num/2)=zeros(Num/2, Num/2);
    H(Num/2+1:Num, 1:Num/2)=zeros(Num/2, Num/2);
    H(Num/2+1:Num, Num/2+1:Num)=zeros(Num/2, Num/2);
    
    fp.H=H;
    fp.rhs=rhs;
    fp.u=u;

elseif strcmp(precond,'gs')
    %reordering
    [fp.Per]=redblack(G);
    
    [MATRICES, fp]=Prec(D, precond, fp, MATRICES);
 
    H=MATRICES.M\MATRICES.N;
    rhs=(fp.Per)'*rhs;
    rhs=(MATRICES.M)\rhs;
    u=(fp.Per)'*u;
    
    fp.H=H;
    fp.M=MATRICES.M;
    fp.N=MATRICES.N;
    fp.rhs=rhs;
    fp.u=u;
    
elseif strcmp(precond,'trisplit')
    
    %reordering
    [fp.Per]=redblack(G);
    [MATRICES, fp]=Prec(D, precond, fp, MATRICES);
    
    H=MATRICES.M\MATRICES.N;
    rhs=(fp.Per)'*rhs;
    rhs=(MATRICES.M)\rhs;
    u=(fp.Per)'*u;
    
    fp.H=H;
    fp.M=MATRICES.M;
    fp.N=MATRICES.N;
    fp.rhs=rhs;
    fp.u=u;

elseif strcmp(precond,'alternating')
    
    [fp.Per]=redblack(G);
    [MATRICES, fp]=Prec(D, 'gs', fp, MATRICES);
    
    H1=MATRICES.M\MATRICES.N;
    rhs1=(fp.Per)'*rhs;
    rhs1=(MATRICES.M)\rhs1;
    
    fp.H1=H1;
    fp.M=MATRICES.M;
    fp.N=MATRICES.N;
    fp.rhs1=rhs1;
    
    [MATRICES, fp]=Prec(D, 'trisplit', fp, MATRICES);
    
    H2=MATRICES.M\MATRICES.N;
    rhs2=(fp.Per)'*rhs;
    rhs2=(MATRICES.M)\rhs2;
    
    fp.H2=H2;
    fp.P=MATRICES.M;
    fp.Q=MATRICES.N;
    fp.rhs2=rhs2;
    
    u=(fp.Per)'*u;
    fp.u=u;
    
    if ( ~ P_regular(fp.M, fp.N) ) || ( ~ P_regular(fp.P, fp.Q) )
        error('P-regularity constraint NOT respected');
    end

    
else
    error('Invalid preconditioner inserted');
end

end