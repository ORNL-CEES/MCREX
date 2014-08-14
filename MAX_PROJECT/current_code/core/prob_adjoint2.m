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
        for i=1:size(P,1)
            for j=1:size(P,2)
                if sum(abs(A(:,i)))>0
                    P(i,j)=(A(j,i)~=0)/(length(find(A(:,i))));
                end
            end
        end

   else
        P=zeros(size(A));
        for i=1:size(P,1)
            for j=1:size(P,2)
                if sum(abs(A(:,i)))>0
                    P(i,j)=abs(A(j,i)).^p2/(sum(abs(A(:,i)).^p2));
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
    
    cdf=P;
    %computation of the cumulative probability
    for i=1:size(cdf,1)
        for j=2:size(cdf,2)
            aux=max(find(cdf(i,1:j-1)));
            if (cdf(i,j)~=0 && ~isempty(aux))
                cdf(i,j)=cdf(i,j)+cdf(i,aux);
            end
        end
    end
    
    P=sparse(P);
    cdf=sparse(cdf);
    
end