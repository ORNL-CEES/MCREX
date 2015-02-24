function [D,iter] = gdd_to_sdd (A, eps, pol)

    iter=0;
    n=size(A,1);
    D=speye((size(A)));
    C=sparse(A);
    absDiag = abs(diag(C));
    
    t=0;
    
    if strcmp(pol, 'row')
        while t<n 
            t=0;
            absDiag = abs(diag(C));
            Sum = sum(abs(C), 2) - absDiag; 
            for i=1:n
                if abs(C(i,i)) > Sum(i)
                    t=t+1;
                end                
            end
            if t<n
                B=zeros(size(C));
                for i=1:n
                    B(i,i) = ( Sum(i)+eps ) / ( abs(C(i,i)) + eps);
                    aux=find( C(:,i) );
                    C(aux,i) = C(aux,i) * B(i,i);
                end
                B=sparse(B);
                D=sparse(D*B);
                iter=iter+1;
            end
        end
    end
    
    if strcmp(pol, 'col')
        while t<n 
            t=0;
            absDiag = abs(diag(C));
            Sum = sum(abs(C), 1)' - absDiag; 
            for i=1:n
                if abs(C(i,i)) > Sum(i)
                    t=t+1;
                end                
            end
            if t<n
                B=zeros(size(C));
                for i=1:n
                    B(i,i) = ( Sum(i)+eps ) / ( abs(C(i,i)) + eps);
                    aux=find( C(i,:) );
                    C(i,aux) = C(i,aux) * B(i,i);
                end
                B=sparse(B);
                D=sparse(B*D);
                iter=iter+1;
            end
        end
    end
end