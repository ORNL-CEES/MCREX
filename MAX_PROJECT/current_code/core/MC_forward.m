function [x, y, X]=MC_forward(A, b, P, cdf, n_walks, max_step)

    X=zeros(n_walks,size(b,1));

   for k=1:size(b,1)
        parfor walk=1:n_walks
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
   end

   %computation of the expected value for the updating vector
   x=mean(X,1)';
   Y=X.^2;
   y=sqrt((mean(Y,1)'-(x.^2))./(size(X,1))); 
end