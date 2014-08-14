function [post, cdf_bay]=bayes(prior, likely)

    display('Bayesian updating');

    if (size(prior)~=size(likely))
        err('prior and likelihood have different dimensions');
    end

    for i=1:size(prior,1)    
        if (find(prior(i,:))~=find(likely(i,:)))
            err('prior and likely matrix are not consistent');
        end
    end

    post=zeros(size(prior));

    for i=1:size(prior,1)
        aux=find(prior(i,:));
        if (~isempty(aux))
            for j=1:size(aux,2)
                post(i,aux(j))=prior(i,aux(j))*likely(i,aux(j))/(sum(prior(i,aux).*likely(i,aux)));
            end
        end
    end

    cdf_bay=post;
    %computation of the cumulative probability
    for i=1:size(cdf_bay,1)
        for j=2:size(cdf_bay,2)
            aux=max(find(cdf_bay(i,1:j-1)))+0;
            if (cdf_bay(i,j)~=0 && ~isempty(aux))
                cdf_bay(i,j)=cdf_bay(i,j)+cdf_bay(i,aux);
            end
        end
    end
    
    post=sparse(post);
    cdf_bay=sparse(cdf_bay);    
    
end