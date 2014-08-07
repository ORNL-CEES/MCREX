function [x, y, TALLY, time]=MC_adjoint_error2(A, b, P, cdf, Pb, cdfb, n_walks, max_step, cutoff)

    L=length(n_walks);
    
    x=[];
    y=[];
    time=[];
    TALLY=[];

    for l=1:L
        
       X=zeros(n_walks(l),size(b,1));
       tally=zeros(size(b))';
       
       start=cputime;
       
       for walk=1:n_walks(l)
           aux=rand;

           %it detects what is the inital status of the chain
           previous=min(find(cdfb>aux));
           W=sum(abs(b))*sign(b(previous));
           Wf=W;
           tally(previous)=tally(previous)+1;
           %W=b(previous)*length(find(P(previous,:)));
           X(walk,previous)=X(walk,previous)+W;
           i=1;
           while i<=max_step && W/Wf>cutoff
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
       finish=cputime;
       %computation of the expected value for the updating vector
       x=[ x mean(X,1)'];
       Y=X.^2;
       y=[ y sqrt((mean(Y,1)'-(mean(X,1)'.^2))./(size(X,1))) ]; 
       time=[time finish-start];
       TALLY=[TALLY; tally];
    end
end
 