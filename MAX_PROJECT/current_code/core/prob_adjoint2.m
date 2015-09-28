function [Pb, cdfb, P, cdf]=prob_adjoint2(A, b, p1, p2)

    display('Building of transition matrix');

    if (p1 == 0)
        Pb=sparse(size(A,1),1);
        for i=1:size(Pb,1)
            Pb(i)=(b(i)~=0)/(length(find(b)));
        end
    else        

        %initial probability distribution
        Pb=sparse(size(A,1),1);
        for i=1:size(A,1)
            Pb(i)=(abs(b(i)).^p1)/(sum(abs(b).^p1));
        end

    end

    
   if (p2 == 0)    
        P=sparse(zeros(size(A)));
        cdf=sparse(zeros(size(A)));
        for i=1:size(A,2)
            Prow_ind=find(A(:,i));
            num_elem=length(Prow_ind);   
            for j=1:1:num_elem
                P(i,Prow_ind(j))=1/num_elem;
                cdf(i,Prow_ind(j))=j/num_elem;
            end       
        end

   else
        P=sparse(zeros(size(A)));
        cdf=sparse(zeros(size(A)));
        for i=1:size(A,2)
            Prow_ind=find(A(:,i));
            num_elem=length(Prow_ind);  
            sump=(sum(abs(A(:,i)).^p2));
            for j=1:1:num_elem
                P(i,Prow_ind(j))=abs(A(Prow_ind(j),i)).^p2/sump;
                cdf(i,Prow_ind(j))=sum(abs(A(Prow_ind(1:j),i)))/sump;
            end        
        end
   end
    
    
    %computation of the cumulative initial probability: 
    cdfb=Pb;
    for i=2:size(cdfb)
        aux=max(find(cdfb(1:i-1)));
        if(cdfb(i)~=0 && ~isempty(aux))
            cdfb(i)=cdfb(i)+cdfb(aux);
        end
    end
    
 
    P=sparse(P);
    cdf=sparse(cdf);
    
end