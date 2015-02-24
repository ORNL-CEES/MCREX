function [Pb, cdfb]=prob_adjoint_rhs(b, p)

display('Updating of the initial probability');

    if p==0
        
        Pb=zeros(size(A,1),1);
        leng=length(b);
        index=find(b);
        n=length(index);
        for i=1:n
            Pb(index(i))=1/leng;
        end
        
    else

        %initial probability distribution
        Pb=zeros(size(b,1),1);
        index=find(b);
        n=length(index);
        b_abs=abs(b).^p;
        sump=sum(b_abs);
        for i=1:n
            Pb(index(i))=b_abs(index(i))/sump;
        end  
        
    end
        
    %computation of the cumulative initial probability: 
    cdfb=Pb;
    for i=2:size(cdfb)
        aux=find(cdfb(1:i-1), 1, 'last' )+0;
        if(cdfb(i)~=0 && ~isempty(aux))
            cdfb(i)=cdfb(i)+cdfb(aux);
        end
    end
end