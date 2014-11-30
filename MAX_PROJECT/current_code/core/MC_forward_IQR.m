function [X, Y, IQR]=MC_forward_IQR(A, b, P, cdf, n_walks, max_step)

    X=zeros(size(b));
    Y=zeros(size(b));
    IQR=cell(size(b,1),1);

   for k=1:size(b,1)
       x=[];
        for walk=1:n_walks
            estim=0;
            previous=k;
            W=1;
            current=k;
            estim=estim+W*b(current);
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
                    estim=estim+W*b(current);
                    i=i+1;
                    previous=current;
                end
                x=[x estim];
                IQR{k}=[IQR{k} iqr(x)/std(x)];
        end
        X(k)=mean(x);
        Y(k)=var(x);  
   end

end