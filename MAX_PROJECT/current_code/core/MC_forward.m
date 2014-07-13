function [x, y, Z]=MC_forward(Z, A, b, P, cdf, n_walks, max_step)

    X=zeros(n_walks,size(b,1));

   for k=1:size(b,1)
        for walk=1:n_walks
        previous=k;
        W=1;
        current=k;
        X(walk,k)=X(walk,k)+W*b(current);
        i=1;
            while i<=max_step
                aux=rand;
                current=min(find(cdf(previous,:)>aux));
                W=W*A(previous,current)/P(previous,current);
                                
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
   Z=[Z; X];
   x=mean(Z,1)';
   Y=Z.^2;
   y=sqrt((mean(Y,1)'-(x.^2))./(size(Z,1))); 
end