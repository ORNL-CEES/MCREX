function [X, VAR, NWALKS]=MC_forward_adapt(A, b, P, cdf, n_walks, max_step, var_cut)

    X=zeros(size(b));
    NWALKS=zeros(size(b));
    VAR=zeros(size(b));

   for k=1:size(b,1)
       count=0;
       x=[];
       ratio=1;
       while count < 1000*size(b,1) && ratio > var_cut
            for walk=1:n_walks
                estim=0;
                previous=k;
                W=1;
                current=k;
                estim=estim+W*b(current);
                i=1;
                while i<=max_step
                    aux=rand;
                    current=min(find(cdf(previous,:)>aux));
                    W=W*A(previous,current)/P(previous,current);

                    if W==0
                        break;
                    end
                    estim=estim+W*b(current);
                    i=i+1;
                    previous=current;
                end
                x=[x; estim];
                dvar=sqrt((mean(x.^2)-(mean(x)).^2)/size(x,1));
                ratio=dvar/abs(mean(x));
            end
            count=count+n_walks;
       end
       X(k)=mean(x);
       NWALKS(k)=count;
       VAR(k)=dvar;
   end

end