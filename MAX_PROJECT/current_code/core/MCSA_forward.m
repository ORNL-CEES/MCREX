function [sol, rel_res, var, VAR, DX, NWALKS, count]=MCSA_forward(fp, P, cdf, numer, stat)

VAR=[];
NWALKS=[];
DX=[];
n_walks=stat.nwalks;
max_step=stat.max_step;
eps=numer.eps;
rich_it=numer.rich_it;

if stat.adapt==1    
    
    var_cut=stat.varcut;
  
    if ~ strcmp(fp.precond, 'alternating')

        if ~ asyn_check(fp.H)
            error('Iteration matrix does not have spectral radius less than 1');
        end 

        display(strcat('Forward MCSA with', {' '},  fp.precond, ' preconditioning started'));

        H=fp.H;
        rhs=fp.rhs;
        sol=ones(size(H,1),1);
        var=zeros(size(H,1),1);

        %matrix to be used for the computation of the redisual at each Richardson
        %iteration
        B=(eye(size(H))-H);

        rel_residual=norm(rhs-B*sol,2)/norm(rhs,2);
        count=1;
        while(rel_residual>eps && count<=rich_it)
            display(strcat('iteration number', {' '}, num2str(count)))
            sol=H*sol+rhs; 
            r=rhs-B*sol;
            [dx, dvar, nwalks]=MC_forward_adapt(H, r, P, cdf, n_walks, max_step, var_cut);
            NWALKS=[NWALKS nwalks];       
            sol=sol+dx;
            r=rhs-B*sol;
            display(strcat('residual norm: ', num2str(norm(r))));    
            rel_residual=norm(r,2)/norm(rhs,2);
            VAR=[VAR dvar];
            DX=[DX dx];
            count=count+1;
        end
        count=count-1;
        rel_res=rel_residual;

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

        rel_residual=norm(rhs1-B1*sol,2)/norm(rhs1,2);
        count=1;
        while(rel_residual>eps && count<=rich_it)
           if mod(count,2)==1
               sol=H1*sol+rhs1; 
               r=rhs1-B1*sol;
               [dx, dvar, nwalks]=MC_forward_adapt(H1, r, P.P1, cdf.cdf1, n_walks, max_step);
               NWALKS=[NWALKS nwalks];
               sol=sol+dx;
               r=rhs1-B1*sol;
               VAR1=[VAR1 dvar];
               VAR=[VAR dvar];
               DX=[DX dx];
               rel_residual=norm(r,2)/norm(rhs1,2);
           else
               sol=H2*sol+rhs2;
               r=rhs2-B2*sol;
               [dx, dvar, nwalks]=MC_forward_adapt(H2, r, P.P2, cdf.cdf2, n_walks, max_step);
               NWALKS=[NWALKS nwalks];
               sol=sol+dx;  
               r=rhs2-B2*sol;
               VAR2=[VAR2 dvar];
               VAR=[VAR dvar];
               DX=[DX dx];
               rel_residual=norm(r,2)/norm(rhs2,2);
           end

           count=count+1;
        end
        count=count-1;
        rel_res=rel_residual;

        for i=1:count
            var=var+abs(H2*((H1*H2)^(count-1))*VAR1(:,i));
            var=var+abs((H1*H2)^(count-1)*VAR2(:,i));
        end

    end
    
else
    
    u=fp.u;
    VAR=[];
    NWALKS=[];
    DX=[];

    if ~ strcmp(fp.precond, 'alternating')

        if ~ asyn_check(fp.H)
            error('Iteration matrix does not have spectral radius less than 1');
        end 

        display(strcat('Forward MCSA with', {' '},  fp.precond, ' preconditioning started'));

        H=fp.H;
        rhs=fp.rhs;
        sol=ones(size(H,1),1);
        var=zeros(size(H,1),1);

        %matrix to be used for the computation of the redisual at each Richardson
        %iteration
        B=(eye(size(H))-H);

        rel_residual=norm(rhs-B*sol,2)/norm(rhs,2);
        count=1;
        while(rel_residual>eps && count<=rich_it)
            display(strcat('iteration number', {' '}, num2str(count)))
            X=[];
            sol=H*sol+rhs; 
            r=rhs-B*sol;
            [dx, dvar, X]=MC_forward(X, H, r, P, cdf, n_walks, max_step);
            while max(dvar)>( norm(dx)/100 )
                fprintf(max(dvar)>( norm(dx)/100 ));
                dnwalks=n_walks;
                [dx, dvar, X]=MC_forward(X, H, r, P, cdf, dnwalks, max_step);
            end
            NWALKS=[NWALKS size(X,1)];
            sol=sol+dx;
            r=rhs-B*sol;
            display(strcat('residual norm: ', num2str(norm(r))));    
            rel_residual=norm(r,2)/norm(rhs,2);
            VAR=[VAR dvar];
            DX=[DX dx];
            count=count+1;
        end
        count=count-1;
        rel_res=rel_residual;

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

        rel_residual=norm(rhs1-B1*sol,2)/norm(rhs1,2);
        count=1;
        while(rel_residual>eps && count<=rich_it)
           X=[];
           if mod(count,2)==1
               sol=H1*sol+rhs1; 
               r=rhs1-B1*sol;
               [dx, dvar, X]=MC_forward(X, H1, r, P.P1, cdf.cdf1, n_walks, max_step);
               while max(dvar)>( norm(dx)/100 )
                   fprintf(max(dvar)>( norm(dx)/100 ));
                   dnwalks=n_walks;
                   [dx, dvar, X]=MC_forward(X, H1, r, P.P1, cdf.cdf1, dnwalks, max_step);
               end
               NWALKS=[NWALKS size(X,1)];
               sol=sol+dx;
               r=rhs1-B1*sol;
               VAR1=[VAR1 dvar];
               VAR=[VAR dvar];
               DX=[DX dx];
               rel_residual=norm(r,2)/norm(rhs1,2);
           else
               sol=H2*sol+rhs2;
               r=rhs2-B2*sol;
               [dx, dvar, X]=MC_forward(X, H2, r, P.P2, cdf.cdf2, n_walks, max_step);
               while max(dvar)>( norm(dx)/100 )
                   dnwalks=n_walks;
                   [dx, dvar, X]=MC_forward(X, H2, r, P.P2, cdf.cdf2, dnwalks, max_step);
               end
               NWALKS=[NWALKS size(X,1)];
               sol=sol+dx;  
               r=rhs2-B2*sol;
               VAR2=[VAR2 dvar];
               VAR=[VAR dvar];
               DX=[DX dx];
               rel_residual=norm(r,2)/norm(rhs2,2);
           end

           count=count+1;
        end
        count=count-1;
        rel_res=rel_residual;

        for i=1:count
            var=var+abs(H2*((H1*H2)^(count-1))*VAR1(:,i));
            var=var+abs((H1*H2)^(count-1)*VAR2(:,i));
        end

    end
    
end

end
