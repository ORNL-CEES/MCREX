function [x, y, err, time]=MC_forward_error2(u_ex, A, b, P, cdf, n_walks, max_step, cutoff)
   
    L=length(n_walks);
    
    err=NaN*ones(L,size(b,1));
    x=[];
    y=[];
    time=[];
    
    for l=1:L
        X=zeros(n_walks(l),size(b,1));
        start=cputime;
       for k=1:size(b,1)
            parfor walk=1:n_walks(l)
            previous=k;
            W=1;
            Wf=W;
            current=k;
            X(walk,k)=X(walk,k)+W*b(current);
            i=1;
                while i<=max_step && W/Wf>cutoff
                    aux=rand;
                    cdfrow_ind=find(cdf(previous,:));
                    if ~isempty(cdfrow_ind)
                        current=min(find(cdf(previous,cdfrow_ind)>aux));
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
            err(l,k)=abs(mean(X(:,k))-u_ex(k))/abs(u_ex(k));
       end
       finish=cputime;
       x=[x mean(X,1)']; 
       Y=X.^2;
       y=[y sqrt((mean(Y,1)'-mean(X,1)'.^2)/(sqrt(n_walks(l))))]; 
       time=[time finish-start];
    end
end