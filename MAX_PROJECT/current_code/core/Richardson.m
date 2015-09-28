function [sol, count, rel_residual]=Richardson(fp, numer)

rich_it=numer.rich_it;
eps=numer.eps;

if ~ strcmp(fp.precond, 'alternating')

%     if ~ asyn_check((fp.H)',P)
%         error('Iteration matrix does not have spectral radius less than 1');
%     end   

    display(strcat('Richardson with', {' '}, fp.precond, ' preconditioning started'));

    H=fp.H;
    rhs=fp.rhs;
    sol=rand(size(H,1),1);

    %matrix to be used for the computation of the residual at each Richardson
    %iteration
    B=sparse((speye(size(H))-H));     

    rel_residual=norm(rhs-B*sol,2)/norm(rhs,2);
    count=0;

    r=rhs-B*sol;
    R=norm(r,2)/norm(rhs,2);
    
    while(rel_residual>eps && count<=rich_it)
             count=count+1;
             display(strcat('iteration number', {' '}, num2str(count)))
             sol=sol+r; 
             r=rhs-B*sol;
             display(strcat('residual norm: ', num2str(norm(r)/norm(rhs))));
             rel_residual=norm(r,2)/norm(rhs,2);
             R=[R rel_residual];          
    end

else 

    if ~ asyn_check(fp.H1)
        error('Iteration matrix H1 does not have spectral radius less than 1');
    end

    if ~ asyn_check(fp.H2)
        error('Iteration matrix H2 does not have spectral radius less than 1');
    end   

    display('Richardson MCSA with alternating method started');

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
    count=0;


    while(rel_residual>eps && count<=rich_it)
        count=count+1;
        if mod(count,2)==1
            sol=sol+r; 
            r=rhs1-B1*sol;
            rel_residual=norm(r,2)/norm(rhs1,2);
        else
            sol=sol+r; 
            r=rhs2-B2*sol;
            rel_residual=norm(r,2)/norm(rhs2,2);
        end
    end

end


