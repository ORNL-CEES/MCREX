function [Pb, cdfb]=prob_adjoint_rhs(b, p)

display('Updating of the inital probability');

    if p==0
        
        Pb=zeros(size(A,1),1);
        for i=1:size(Pb,1)
            Pb(i)=(b(i)~=0)/(length(find(b)));
        end
        
    else

        %initial probability distribution
        Pb=zeros(size(b,1),1);
        for i=1:size(Pb,1)
            Pb(i)=(abs(b(i)).^p)/(sum(abs(b)).^p);
        end  
        
    end
        
    %computation of the cumulative initial probability: 
    cdfb=Pb;
    for i=2:size(cdfb)
        aux=max(find(cdfb(1:i-1)))+0;
        if(cdfb(i)~=0 && ~isempty(aux))
            cdfb(i)=cdfb(i)+cdfb(aux);
        end
    end

end