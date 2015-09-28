function [sol, rel_res, VAR, REL_RES, DX, NWALKS, tally, count, reject]=SEQ_adjoint2(fp, dist, P, cdf, numer, stat)

VAR=[];
DX=[];
NWALKS=[];
tally=[];
REL_RES=[];
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
    sol=rand(size(H,1),1);

    %matrix to be used for the computation of the redisual at each Richardson
    %iteration
    B=sparse(speye(size(H))-H);

    r=rhs-B*sol;
    rel_residual=norm(rhs-B*sol,2)/norm(rhs,2);
    REL_RES=[REL_RES rel_residual];
    count=1;

    if stat.adapt_walks==1 || stat.adapt_cutoff==1
        while(rel_residual>eps && count<=rich_it)
            display(strcat('iteration number', {' '}, num2str(count)))
            n_walks=length(find(r));
            [Pb, cdfb]=prob_adjoint_rhs(r, dist);
            
             if stat.adapt_walks==1 && count>1 && REL_RES(end)>REL_RES(end-1)
                 stat.varcut=stat.varcut/ceil(REL_RES(end)/min(REL_RES));
             end                
            
            [dx, dvar, tot_walks, resid, aux, rej]=MC_adjoint_adapt2(H, r, P, cdf, Pb, cdfb, stat);
            NWALKS=[NWALKS tot_walks];
            tally=[tally aux];
            reject=[reject rej];
            sol=sol+dx;
            r=rhs-B*sol;
            display(strcat('residual norm: ', num2str(norm(r,2)/norm(rhs,2))));        
            rel_residual=norm(r,2)/norm(rhs,2);
            REL_RES=[REL_RES rel_residual];
            VAR{count}=dvar;
            RES{count}=resid;
            DX=[DX dx];
            count=count+1;
        end
    else
        while(rel_residual>eps && count<=rich_it)
            display(strcat('iteration number', {' '}, num2str(count)))
            [Pb, cdfb]=prob_adjoint_rhs(r, dist);
            
             if stat.adapt_walks==1 && count>1 && REL_RES(end)>REL_RES(end-1)
                 stat.varcut=stat.varcut/ceil(REL_RES(end)/min(REL_RES));
             end                
            
            [dx, dvar, tot_walks, ~]=MC_adjoint(H, r, P, cdf, Pb, cdfb, n_walks, max_step);
            NWALKS=[NWALKS tot_walks];
            sol=sol+dx;
            r=rhs-B*sol;
            display(strcat('residual norm: ', num2str(norm(r,2)/norm(rhs,2))));        
            rel_residual=norm(r,2)/norm(rhs,2);
            REL_RES=[REL_RES rel_residual];
            VAR=[VAR dvar];
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

    display('Adjoint SEQ with alternating method started');

    H1=fp.H1;
    H2=fp.H2;
    P1=P.P1;
    P2=P.P2;
    cdf1=cdf.cdf1;
    cdf2=cdf.cdf2;
    B1=sparse(speye(size(H1))-H1);
    B2=sparse(speye(size(H2))-H2);

    rhs1=fp.rhs1;
    rhs2=fp.rhs2;
    sol=rand(size(H,1),1);
   
    r=rhs1-B1*sol;
    rel_residual=norm(rhs1-B1*sol,2)/norm(rhs1,2);
    REL_RES=[REL_RES rel_residual];
    count=1;

    if stat.adapt_walks==1 || stat.adapt_cutoff==1
        while(rel_residual>eps && count<=rich_it)
            if mod(count,2)==1
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                
                if stat.adapt_walks==1 && count>1 && REL_RES(end)>REL_RES(end-1)
                   stat.varcut=stat.varcut/ceil(REL_RES(end)/min(REL_RES));
                end                    
                
                [dx, dvar, resid, tot_walks, aux, rej]=MC_adjoint_adapt2(H1, r, P.P1, cdf.cdf1, Pb, cdfb, stat);
                NWALKS=[NWALKS tot_walks];
                tally=[tally aux];
                reject=[reject rej];                
                sol=sol+dx;
                r=rhs1-B1*sol;
                VAR{count}=dvar;
                RES{count}=resid;
                DX=[DX dx]; 
                rel_residual=norm(r,2)/norm(rhs1,2);
                REL_RES=[REL_RES rel_residual];
            else
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                
                if stat.adapt_walks==1 && count>1 && REL_RES(end)>REL_RES(end-1)
                   stat.varcut=stat.varcut/ceil(REL_RES(end)/min(REL_RES));
                end                 
                
                [dx, dvar, resid, tot_walks, aux, rej]=MC_adjoint_adapt2(H2, r, P.P2, cdf.cdf2, Pb, cdfb, stat);
                NWALKS=[NWALKS tot_walks];
                tally=[tally aux];
                reject=[reject rej];           
                sol=sol+dx; 
                r=rhs2-B2*sol; 
                VAR{count}=dvar;
                RES{count}=resid;
                DX=[DX dx];
                rel_residual=norm(r,2)/norm(rhs2,2);
                REL_RES=[REL_RES rel_residual];
            end

            count=count+1;
        end
    else
        while(rel_residual>eps && count<=rich_it)
            X=[];
            if mod(count,2)==1
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);  
                [dx, dvar, tot_walks, ~]=MC_adjoint(H1, r, P1, cdf1, Pb, cdfb, n_walks, max_step);
                NWALKS=[NWALKS tot_walks];
                sol=sol+dx;
                r=rhs1-B1*sol;
                VAR=[VAR dvar];
                DX=[DX dx];
                rel_residual=norm(r,2)/norm(rhs1,2);
                REL_RES=[REL_RES rel_residual];
            else
                [Pb, cdfb]=prob_adjoint_rhs(r, dist);
                [dx, dvar, tot_walks, ~]=MC_adjoint(H2, r, P2, cdf2, Pb, cdfb, n_walks, max_step);
                NWALKS=[NWALKS tot_walks];
                sol=sol+dx; 
                r=rhs2-B2*sol; 
                VAR=[VAR dvar];
                DX=[DX dx];
                rel_residual=norm(r,2)/norm(rhs2,2);
                REL_RES=[REL_RES rel_residual];
            end

            count=count+1;
        end
    end

    count=count-1;
    rel_res=rel_residual; 

end

end