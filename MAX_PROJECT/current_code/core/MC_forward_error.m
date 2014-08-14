function [x, y, err]=MC_forward_error(u_ex, A, b, P, cdf, n_walks, max_step)
   
    L=length(n_walks);
    
    err=NaN*ones(L,size(b,1));
    x=[];
    y=[];
    
    for l=1:L
        X=zeros(n_walks(l),size(b,1));
       for k=1:size(b,1)
            for walk=1:n_walks(l)
            previous=k;
            W=1;
            current=k;
            X(walk,k)=X(walk,k)+W*b(current);
            i=1;
                while i<=max_step
                    aux=rand;
                    if sum(abs(P(previous,:)))>0                    
                        current=min(find(cdf(previous,:)>aux));
                        W=W*A(previous,current)/P(previous,current);
                    else
                        W=0;
                    end
                    if W==0
                        break;
                    end
                    X(walk,k)=X(walk,k)+W*b(current);
                    i=i+1;
                    previous=current;
                end
            end
            err(l,k)=abs(mean(X(:,k))-u_ex(k))/abs(u_ex(k));
       end
       x=[x mean(X,1)']; 
       Y=X.^2;
       y=[y sqrt((mean(Y,1)'-mean(X,1)'.^2)/(sqrt(n_walks(l))))]; 
    end
end