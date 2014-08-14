function [P, cdf]=prob_forward(A, p)

display('Building of transition matrix');

    if (p == 0)
        P=zeros(size(A));
        for i=1:size(P,1)
            if isempty(A(i,:))
                break;
            else            
                for j=1:size(P,2)
                    P(i,j)=(A(i,j)~=0)/(length(find(A(i,:))));
                end
            end
        end

    else

        P=zeros(size(A));
        for i=1:size(P,1)
            if isempty(A(i,:))
                break;
            else            
                for j=1:size(P,2)
                    P(i,j)=abs(A(i,j)).^p/(sum(abs(A(i,:)).^p));
                end
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
