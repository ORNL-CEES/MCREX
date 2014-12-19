function [Pb, cdfb, P, cdf]=prob_adjoint2(A, b, p1, p2)

    display('Building of transition matrix');

    if (p1 == 0)
        Pb=zeros(size(A,1),1);
        for i=1:size(Pb,1)
            Pb(i)=(b(i)~=0)/(length(find(b)));
        end
    else        

        %initial probability distribution
        Pb=zeros(size(A,1),1);
        for i=1:size(A,1)
            Pb(i)=(abs(b(i)).^p1)/(sum(abs(b).^p1));
        end

    end

    
   if (p2 == 0)    
        P=zeros(size(A));
        cdf=zeros(size(A));
        for i=1:size(P,1)
            Prow_ind=find(A(:,i));
            num_elem=length(find(A(:,i)));   
            if ~isempty(find(A(:,i)))
               aux=find(A(:,i));
                for j=1:1:num_elem
                    P(i,Prow_ind(j))=(A(Prow_ind(j),i)~=0)/num_elem;
                    cdf(i,Prow_ind(j))=j/num_elem;
                end
            end
        end

   else
        P=zeros(size(A));
        cdf=zeros(size(A));
        for i=1:size(P,1)
           Prow_ind=find(A(:,i));
            num_elem=length(find(A(:,i)));  
            sump=(sum(abs(A(:,i)).^p));
            for j=1:1:num_elem
                if ~isempty(find(A(:,i)))
                    P(i,Prow_ind(j))=abs(A(Prow_ind(j),i)).^p/sump;
                    cdf(i,Prow_ind(j))=sum(abs(A(Prow_ind(1:j),i)))/sump;
                end
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