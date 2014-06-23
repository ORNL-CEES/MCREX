function [sol, rel_err, var, VAR, count]=MCSA_forward(fp, P, cdf, rich_it, n_walks, max_step, eps)

u=fp.u;
VAR=[];

if ~ strcmp(fp.precond, 'alternating')
    
    if ~ asyn_check(fp.H)
        error('Iteration matrix does not have spectral radius less than 1');
    end 
    
    display(strcat('Forward MCSA with ', fp.precond, ' preconditioning started'));

    H=fp.H;
    rhs=fp.rhs;
    sol=ones(size(H,1),1);
    var=zeros(size(H,1),1);

    %matrix to be used for the computation of the redisual at each Richardson
    %iteration
    B=(eye(size(H))-H);

    rel_error=norm(u-sol,2)/norm(u,2);
    count=1;
    while(rel_error>eps && count<=rich_it)
       sol=H*sol+rhs; 
       r=rhs-B*sol;
       [dx, dvar]=MC_forward(H, r, P, cdf, n_walks, max_step);
       sol=sol+dx;
       rel_error=norm(u-sol,2)/norm(u,2);
       VAR=[VAR dvar];
       count=count+1;
    end
    count=count-1;
    rel_err=rel_error;
    
    for i=1:count
        var=var+abs(H^(count-i)*VAR(:,i));
    end
   
else
    
    if ~ asyn_check(fp.H1)
        error('Iteration matrix H1 does not have spectral radius less than 1');
    end

    if ~ asyn_check(fp.H2)
        error('Iteration matrix H2 does not have spectral radius less than 1');
    end
    
    display('Forward MCSA with alternating method started');
    
    H1=fp.H1;
    H2=fp.H2;
    rhs1=fp.rhs1;
    rhs2=fp.rhs2;
     
    sol=ones(size(H1,1),1);
    var=zeros(size(H1,1),1);
    
    VAR1=[];
    VAR2=[];
    
    %matrix to be used for the computation of the redisual at each Richardson
    %iteration
    B1=(eye(size(H1))-H1);
    B2=(eye(size(H2))-H2);    
    
    rel_error=norm(u-sol,2)/norm(u,2);
    count=1;
    while(rel_error>eps && count<=rich_it)
       sol=H1*sol+rhs1; 
       r=rhs1-B1*sol;
       [dx, dvar]=MC_forward(H1, r, P.P1, cdf.cdf1, n_walks, max_step);
       sol=sol+dx;
       VAR1=[VAR1 dvar];
       sol=H2*sol+rhs2;
       r=rhs2-B2*sol;
       [dx, dvar]=MC_forward(H2, r, P.P2, cdf.cdf2, n_walks, max_step);
       sol=sol+dx;  
       VAR2=[VAR2 dvar];
       rel_error=norm(u-sol,2)/norm(u,2);
       count=count+1;
    end
    count=count-1;
    rel_err=rel_error;
    
    for i=1:count
        var=var+abs(H2*((H1*H2)^(count-1))*VAR1(:,i));
        var=var+abs((H1*H2)^(count-1)*VAR2(:,i));
    end

end