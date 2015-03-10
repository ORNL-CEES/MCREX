function [Pb, cdfb, P, cdf]=prob_adjoint(A, b, p)

display('Building of transition matrix');

    if (p == 0)
        
        Pb=zeros(size(A,1),1);
        for i=1:size(Pb,1)
            Pb(i)=(b(i)~=0)/(length(find(b)));
        end
       
        P=zeros(size(A));
        cdf=zeros(size(A));
        for i=1:size(P,1)
            Prow_ind=find(A(:,i));
            num_elem=length(Prow_ind);            
            for j=1:1:num_elem
                P(i,Prow_ind(j))=(A(Prow_ind(j),i)~=0)/num_elem;
                cdf(i,Prow_ind(j))=j/num_elem;
            end
        end
        
    else

        %initial probability distribution
        Pb=zeros(size(A,1),1);
        for i=1:size(A,1)
            Pb(i)=(abs(b(i)).^p)/(sum(abs(b)).^p);
        end


        P=zeros(size(A));
        cdf=zeros(size(A));
        for i=1:size(P,1)
            Prow_ind=find(A(:,i));
            num_elem=length(find(A(:,i)));  
            sump=(sum(abs(A(:,i)).^p));
            for j=1:1:num_elem
                P(i,Prow_ind(j))=abs(A(Prow_ind(j),i)).^p/sump;
                cdf(i,Prow_ind(j))=sum(abs(A(Prow_ind(1:j),i)))/sump;
            end
        end
    
    end
    
    %computation of the cumulative initial probability: 
    cdfb=Pb;
    for i=2:size(cdfb)
        aux=find(cdfb(1:i-1), 1, 'last' );
        if(cdfb(i)~=0 && ~isempty(aux))
            cdfb(i)=cdfb(i)+cdfb(aux);
        end
    end
    
    
    P=sparse(P);
    cdf=sparse(cdf);
    
end