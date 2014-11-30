function [sol, rel_res, IQR, RES, DX, NWALKS, tally, count]=MCSA_adjoint_IQR(fp, dist, P, cdf, numer, stat)

DX=[];
NWALKS=[];
tally=[];
n_walks=stat.nwalks;
max_step=stat.max_step;
eps=numer.eps;
rich_it=numer.rich_it;
IQR=cell(rich_it,1);
RES=cell(rich_it,1);

if ~ strcmp(fp.precond, 'alternating')

    if ~ asyn_check(fp.H)
        error('Iteration matrix does not have spectral radius less than 1');
    end   

    display(strcat('Adjoint MCSA with', {' '}, fp.precond, ' preconditioning started'));

    H=fp.H;
    rhs=fp.rhs;
    sol=ones(size(H,1),1);

    %matrix to be used for the computation of the redisual at each Richardson
    %iteration
    B=(eye(size(H))-H);    

    rel_residual=norm(rhs-B*sol,2)/norm(rhs,2);
    count=1;
  
    r=rhs-B*sol;
    R=norm(r,2)/norm(rhs,2);
    
    if stat.adapt_walks==1 || stat.adapt_cutoff==1
        while(rel_residual>eps && count<=rich_it)
             display(strcat('iteration number', {' '}, num2str(count)))
             sol=sol+r; 
             r=rhs-B*sol;
             [Pb, cdfb]=prob_adjoint_rhs(r, dist);
             
             if count>1 && norm(R(:,end))/norm(rhs) > norm(R(:,end-1))/norm(rhs)
                 stat.iqr=stat.iqr/10;
             end
                          
             [dx, Iqr, X, aux, resid]=MC_adjoint_adapt_IQR(H, r, P, cdf, Pb, cdfb, stat);
             NWALKS=[NWALKS size(X,1)];
             tally=[tally aux];
             sol=sol+dx;
             r=rhs-B*sol;
             R=[R r];
             display(strcat('residual norm: ', num2str(norm(r)/norm(rhs))));
             rel_residual=norm(r,2)/norm(rhs,2);
             IQR{count}=Iqr;
             RES{count}=resid;
             DX=[DX dx];
             count=count+1;
        end
    else
         while(rel_residual>eps && count<=rich_it)
             display(strcat('iteration number', {' '}, num2str(count)))
             sol=sol+r; 
             r=rhs-B*sol;
             [Pb, cdfb]=prob_adjoint_rhs(r, dist);
             [dx, dvar, X, aux]=MC_adjoint(H, r, P, cdf, Pb, cdfb, n_walks, max_step);
             NWALKS=[NWALKS size(X,1)];
             tally=[tally aux];
             sol=sol+dx;
             r=rhs-B*sol;
             display(strcat('residual norm: ', num2str(norm(r)/norm(rhs))));
             rel_residual=norm(r,2)/norm(rhs,2);         
             DX=[DX dx];
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

    display('Adjoint MCSA with alternating method started');

    H1=fp.H1;
    H2=fp.H2;
    rhs1=fp.rhs1;
    rhs2=fp.rhs2;
    sol=ones(size(H1,1),1);

    %matrix to be used for the computation of the redisual at each Richardson
    %iteration
    B1=(eye(size(H1))-H1);
    B2=(eye(size(H2))-H2);

    rel_residual=norm(rhs1-B1*sol,2)/norm(rhs1,2);
    R=rel_residual;
    
    count=1;
    
    IQR1=[];
    IQR2=[];


    if stat.adapt_walks==1 || stat.adapt_cutoff==1
        while(rel_residual>eps && count<=rich_it)
            if mod(count,2)==1
                sol=sol+r; 
                r=rhs1-B1*sol;
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                
                 if count>1 && norm(R(:,end))/norm(rhs1) > norm(R(:,end-1))/norm(rhs1)
                     stat.iqr=stat.iqr/10;
                 end                               
                
                [dx, Iqr, X, aux, resid]=MC_adjoint_adapt_IQR(H1, r, P.P1, cdf.cdf1, Pb, cdfb, stat);
                NWALKS=[NWALKS size(X,1)];
                tally=[tally aux];
                sol=sol+dx;
                r=rhs1-B1*sol;
                IQR{count}=Iqr;
                RES{count}=resid;
                DX=[DX dx];
                rel_residual=norm(r,2)/norm(rhs1,2);
                R=[R rel_residual];
            else
                sol=sol+r; 
                r=rhs2-B2*sol;
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                
                 if count>1 && norm(R(:,end))/norm(rhs1) > norm(R(:,end-1))/norm(rhs1)
                     stat.iqr=stat.iqr/10;
                 end                  
                
                [dx, Iqr, X, aux, resid]=MC_adjoint_adapt_IQR(H2, r, P.P2, cdf.cdf2, Pb, cdfb, stat);
                NWALKS=[NWALKS size(X,1)];
                tally=[tally aux];
                sol=sol+dx;
                r=rhs2-B2*sol;
                IQR{count}=Iqr;
                RES{count}=resid;
                DX=[DX dx];
                rel_residual=norm(r,2)/norm(rhs2,2);
                R=[R rel_residual];
            end
            count=count+1;
        end
        
    else
        while(rel_residual>eps && count<=rich_it)
            X=[];
            if mod(count,2)==1
                sol=sol+r; 
                r=rhs1-B1*sol;
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                [dx, dvar, X, aux]=MC_adjoint(H1, r, P.P1, cdf.cdf1, Pb, cdfb, n_walks, max_step);
                NWALKS=[NWALKS size(X,1)];
                tally=[tally aux];
                sol=sol+dx;
                r=rhs1-B1*sol;
                DX=[DX dx];
                rel_residual=norm(r,2)/norm(rhs1,2);
            else
                sol=sol+r; 
                r=rhs2-B2*sol;
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                [dx, dvar, X, aux]=MC_adjoint(H2, r, P.P2, cdf.cdf2, Pb, cdfb, n_walks, max_step);
                NWALKS=[NWALKS size(X,1)];
                tally=[tally aux];
                sol=sol+dx;
                r=rhs2-B2*sol;
                DX=[DX dx];
                rel_residual=norm(r,2)/norm(rhs2,2);
            end
            count=count+1;
        end

    end
               
    count=count-1;
    rel_res=rel_residual;

end

end


