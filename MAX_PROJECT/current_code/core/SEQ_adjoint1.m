function [sol, rel_res, var, VAR, DX, NWALKS, tally, count]=SEQ_adjoint1(fp, dist, P, cdf, numer, stat)

VAR=[];
DX=[];
NWALKS=[];
tally=[];
n_walks=stat.nwalks;
max_step=stat.max_step;
eps=numer.eps;
rich_it=numer.rich_it;
  

if ~ strcmp(fp.precond, 'alternating')

    if ~ asyn_check(fp.H)
        error('Iteration matrix does not have spectral radius less than 1');
    end  

    display(strcat('Adjoint SEQ with ', {' '}, fp.precond, ' preconditioning started'));    

    H=fp.H;
    rhs=fp.rhs;
    sol=ones(size(H,1),1);
    var=zeros(size(H,1),1);

    %matrix to be used for the computation of the redisual at each Richardson
    %iteration
    B=(eye(size(H))-H);

    rel_residual=norm(rhs-B*sol,2)/norm(rhs,2);
    count=1;

    if stat.adapt_walks==1 || stat.adapt_cutoff==1
        while(rel_residual>eps && count<=rich_it)
            display(strcat('iteration number', {' '}, num2str(count)))
            X=[];
            r=rhs-B*sol;
            n_walks=length(find(r));
            [Pb, cdfb]=prob_adjoint_rhs(r, dist);
            [dx, dvar, X, aux]=MC_adjoint_adapt(H, r, P, cdf, Pb, cdfb, stat);
            NWALKS=[NWALKS size(X,1)];
            tally=[tally aux];
            sol=sol+dx;
            r=rhs-B*sol;
            display(strcat('residual norm: ', num2str(norm(r,2)/norm(rhs,2))));        
            rel_residual=norm(r,2)/norm(rhs,2);
            VAR=[VAR dvar];
            DX=[DX dx];
            var=var+abs(dvar);
            count=count+1;
        end
    else
        while(rel_residual>eps && count<=rich_it)
            display(strcat('iteration number', {' '}, num2str(count)))
            X=[];
            r=rhs-B*sol;
            [Pb, cdfb]=prob_adjoint_rhs(r, dist);
            [dx, dvar, X, ~]=MC_adjoint(H, r, P, cdf, Pb, cdfb, n_walks, max_step);
            NWALKS=[NWALKS size(X,1)];
            sol=sol+dx;
            r=rhs-B*sol;
            display(strcat('residual norm: ', num2str(norm(r,2)/norm(rhs,2))));        
            rel_residual=norm(r,2)/norm(rhs,2);
            VAR=[VAR dvar];
            DX=[DX dx];
            var=var+abs(dvar);
            count=count+1;
        end
    end


    count=count-1;
    rel_res=rel_residual; 

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
    rel_residual=norm(rhs1-B1*sol,2)/norm(rhs1,2);
    count=1;

    if stat.adapt_walks==1 || stat.adapt_cutoff==1
        while(rel_residual>eps && count<=rich_it)
            X=[];
            if mod(count,2)==1
                r=rhs1-B1*sol;
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                [dx, dvar, X, aux]=MC_adjoint_adapt(H1, r, P1, cdf1, Pb, cdfb, stat);
                NWALKS=[NWALKS size(X,1)];
                tally=[tally aux];
                sol=sol+dx;
                r=rhs1-B1*sol;
                VAR=[VAR dvar];
                DX=[DX dx];
                var=var+abs(dvar);
                rel_residual=norm(r,2)/norm(rhs1,2);
            else
                r=rhs2-B2*sol; 
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                [dx, dvar, X, aux]=MC_adjoint_adapt(H2, r, P2, cdf2, Pb, cdfb, stat);
                NWALKS=[NWALKS size(X,1)];
                tally=[tally aux];
                sol=sol+dx; 
                r=rhs2-B2*sol; 
                VAR=[VAR dvar];
                DX=[DX dx];
                var=var+abs(dvar);
                rel_residual=norm(r,2)/norm(rhs2,2);
            end

            count=count+1;
        end
    else
        while(rel_residual>eps && count<=rich_it)
            X=[];
            if mod(count,2)==1
                r=rhs1-B1*sol;
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                [dx, dvar, X, ~]=MC_adjoint(H1, r, P1, cdf1, Pb, cdfb, n_walks, max_step);
                NWALKS=[NWALKS size(X,1)];
                sol=sol+dx;
                r=rhs1-B1*sol;
                VAR=[VAR dvar];
                DX=[DX dx];
                var=var+abs(dvar);
                rel_residual=norm(r,2)/norm(rhs1,2);
            else
                r=rhs2-B2*sol; 
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                [dx, dvar, X, ~]=MC_adjoint(H2, r, P2, cdf2, Pb, cdfb, n_walks, max_step);
                NWALKS=[NWALKS size(X,1)];
                sol=sol+dx; 
                r=rhs2-B2*sol; 
                VAR=[VAR dvar];
                DX=[DX dx];
                var=var+abs(dvar);
                rel_residual=norm(r,2)/norm(rhs2,2);
            end

            count=count+1;
        end
    end

    count=count-1;
    rel_res=rel_residual; 

end

end