function [x, y]=MC_adjoint(A,b,P,cdf, Pb, cdfb, n_walks, max_step)

    X=zeros(n_walks,size(b,1));

   for walk=1:n_walks
    aux=rand;
    
    %it detects what is the inital status of the chain
    previous=min(find(cdfb>aux));
    W=sum(abs(b))*sign(b(previous));
    X(walk,previous)=X(walk,previous)+W;
    i=1;
        while i<=max_step
              aux=rand;
              if sum(abs(P(previous,:)))==1
                current=min(find(cdf(previous,:)>aux));
                W=W*A(current,previous)/P(previous,current);
              else
                 W=0;
              end
              if W==0
                  % if I reach a null state, I will never move away from
                  % there and the related permutation will not give any
                  % contribution to the estimation of the updating vector,
                  % thus this permutation can be neglected
                 break;
              end
              
              % I insert a contribution on the estimator related to the
              % current state
              X(walk,current)=X(walk,current)+W;
              i=i+1;
              previous=current;
        end
    end
    
    %computation of the expected value for the updating vector
    x=mean(X,1)';
    Y=X.^2;
    y=sqrt((mean(Y,1)'-x.^2)./(n_walks)); 
end
 