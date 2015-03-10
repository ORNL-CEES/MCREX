function [x, y, tot_walks]=MC_forward(A, b, P, cdf, n_walks, max_step)

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
                cdfrow_ind=find(cdf(previous,:));
                if ~isempty(cdfrow_ind)
                    current=find(cdf(previous,cdfrow_ind)>aux, 1 );
                    current=cdfrow_ind(current);
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
    tot_walks=n_walks;
   %computation of the expected value for the updating vector
   x=mean(X,1)';
   y=sqrt(var(X,1)'./(size(X,1))); 
end