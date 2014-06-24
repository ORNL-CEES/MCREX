function [sol, rel_err, var, NWALKS, count]=SEQ_adjoint(fp, dist, P, cdf, rich_it, n_walks, max_step, eps)

u_ex=fp.u;
VAR=[];
NWALKS=[];

if ~ strcmp(fp.precond, 'alternating')
    
    if ~ asyn_check(fp.H)
        error('Iteration matrix does not have spectral radius less than 1');
    end  
    
    display(strcat('Adjoint SEQ with ',fp.precond, ' preconditioning started'));    
    
    H=fp.H;
    rhs=fp.rhs;
    sol=ones(size(H,1),1);
    var=zeros(size(H,1),1);
    rel_error=norm(u_ex-sol,2)/norm(u_ex,2);
    count=1;

    %matrix to be used for the computation of the redisual at each Richardson
    %iteration
    B=(eye(size(H))-H);

    while(rel_error>eps && count<=rich_it)
        X=[];
        r=rhs-B*sol;
        [Pb, cdfb]=prob_adjoint_rhs(r, dist);
        [dx, dvar, X]=MC_adjoint(X, H, r, P, cdf, Pb, cdfb, n_walks, max_step);
        while max(dvar)>( norm(dx)/100 )
            dnwalks=n_walks;
            [dx, dvar, X]=MC_adjoint(X, H, r, P, cdf, Pb, cdfb, dnwalks, max_step);
        end
        NWALKS=[NWALKS size(X,1)];
        sol=sol+dx;
        rel_error=norm(u_ex-sol,2)/norm(u_ex,2);
        VAR=[VAR dvar];
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
    
    display('Adjoint SEQ with alternating method started');
    
    H1=fp.H1;
    H2=fp.H2;
    P1=P.P1;
    P2=P.P2;
    cdf1=cdf.cdf1;
    cdf2=cdf.cdf2;
    B1=(eye(size(H1))-H1);
    B2=(eye(size(H2))-H2);
    
    rhs1=fp.rhs1;
    rhs2=fp.rhs2;
    sol=ones(size(H1,1),1);
    var=zeros(size(H1,1),1);
    rel_error=norm(u_ex-sol,2)/norm(u_ex,2);
    count=1;

    while(rel_error>eps && count<=rich_it)
        X=[];
        if mod(count,2)==1
            r=rhs1-B1*sol;
            [Pb, cdfb]=prob_adjoint_rhs(r, dist);
            [dx, dvar]=MC_adjoint(H1, r, P1, cdf1, Pb, cdfb, n_walks, max_step);
            while max(dvar)>( norm(dx)/100 )
                dnwalks=n_walks;
                [dx, dvar, X]=MC_adjoint(X, H1, r, P1, cdf1, Pb, cdfb, dnwalks, max_step);
            end
            NWALKS=[NWALKS size(X,1)];
            sol=sol+dx;
            var=var+abs(dvar);
            rel_error=norm(u_ex-sol,2)/norm(u_ex,2);
        else
            r=rhs2-B2*sol; 
            [Pb, cdfb]=prob_adjoint_rhs(r, dist);
            [dx, dvar]=MC_adjoint(H1, r, P2, cdf2, Pb, cdfb, n_walks, max_step);
            while max(dvar)>( norm(dx)/100 )
                dnwalks=n_walks;
                [dx, dvar, X]=MC_adjoint(X, H2, r, P2, cdf2, Pb, cdfb, dnwalks, max_step);
            end  
            NWALKS=[NWALKS size(X,1)];
            sol=sol+dx;  
            var=var+abs(dvar);
            rel_error=norm(u_ex-sol,2)/norm(u_ex,2);
        end
        
        count=count+1;
    end
 
    count=count-1;
    rel_err=rel_error; 

end