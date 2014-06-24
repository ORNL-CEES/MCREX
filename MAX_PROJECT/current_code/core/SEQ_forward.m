function [sol, rel_err, var, NWALKS, count]=SEQ_forward(fp, P, cdf, rich_it, n_walks, max_step, eps)

VAR=[];
NWALKS=[];

if ~ strcmp(fp.precond, 'alternating')
    
    if ~ asyn_check(fp.H)
        error('Iteration matrix does not have spectral radius less than 1');
    end
    
    display(strcat('Forward SEQ with ', fp.precond, ' preconditioning started'));
    
    u=fp.u;
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
       X=[];
       r=rhs-B*sol;
       [dx, dvar, X]=MC_forward(X, H, r, P, cdf, n_walks, max_step);
       while max(dvar)>( norm(dx)/100 )
           dnwalks=n_walks;
           [dx, dvar, X]=MC_adjoint(X, H, r, P, cdf, Pb, cdfb, dnwalks, max_step);
       end
       NWALKS=[NWALKS size(X,1)];       
       sol=sol+dx;
       rel_error=norm(u-sol,2)/norm(u,2);
       VAR=[VAR abs(dvar)];
       var=var+abs(dvar);
       count=count+1;
    end
    count=count-1;
    rel_err=rel_error;    
        
else
    
    if ~ asyn_check(fp.H1)
        error('Iteration matrix H1 does not have spectral radius less than 1');
    end

    if ~ asyn_check(fp.H2)
        error('Iteration matrix H2 does not have spectral radius less than 1');
    end
    
    display('Forward SEQ with alternating method started');
    
    u=fp.u;
    H1=fp.H1;
    H2=fp.H2;
    B1=(eye(size(H1))-H1);
    B2=(eye(size(H2))-H2);
    
    rhs1=fp.rhs1;
    rhs2=fp.rhs2;
    sol=ones(size(H1,1),1);
    rel_error=sqrt(sum((u-sol).^2))/sqrt(sum((u).^2));
    count=1;
    
    var=zeros(size(H1,1),1);

    while(rel_error>eps && count<=rich_it) 
        X=[];
        if mod(count,2)==1
            r=rhs1-B1*sol;
            [dx, dvar, X]=MC_forward(X, H1, r, P.P1, cdf.cdf1, n_walks, max_step);
            while max(dvar)>( norm(dx)/100 )
                dnwalks=n_walks;
                [dx, dvar, X]=MC_forward(X, H1, r, P, cdf, Pb, cdfb, dnwalks, max_step);
            end
            NWALKS=[NWALKS size(X,1)];        
            sol=sol+dx;
            var=var+abs(dvar);
            rel_error=norm(u-sol,2)/norm(u,2);
        else
            r=rhs2-B2*sol;
            [dx, dvar, X]=MC_forward(X, H2, r, P.P2, cdf.cdf2, n_walks, max_step);
            while max(dvar)>( norm(dx)/100 )
                dnwalks=n_walks;
                [dx, dvar, X]=MC_forward(X, H2, r, P.P2, cdf.cdf2, Pb, cdfb, dnwalks, max_step);
            end
            NWALKS=[NWALKS size(X,1)];
            sol=sol+dx;
            var=var+abs(dvar);
            rel_error=norm(u-sol,2)/norm(u,2);
        end
        
        count=count+1;
    end
    count=count-1;
    rel_err=rel_error;    
        
end

