function [P, cdf]=prob_forward(A, dist)

    if strcmp(dist, 'UM')
        P=zeros(size(A));
        for i=1:size(P,1)
            for j=1:size(P,2)
                P(i,j)=(A(i,j)~=0)/(length(find(A(i,:))));
            end
        end

               
    elseif strcmp(dist, 'MAO')
        P=zeros(size(A));
        for i=1:size(P,1)
            for j=1:size(P,2)
                P(i,j)=abs(A(i,j))/(sum(abs(A(i,:))));
            end
        end
    end

    cdf=P;
    %computation of the cumulative probability
    for i=1:size(cdf,1)
        for j=2:size(cdf,2)
            aux=max(find(cdf(i,1:j-1)))+0;
            if (cdf(i,j)~=0 && ~isempty(aux))
                cdf(i,j)=cdf(i,j)+cdf(i,aux);
            end
        end
    end

end
