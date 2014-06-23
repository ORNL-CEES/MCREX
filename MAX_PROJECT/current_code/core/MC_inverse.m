function [inv_A]=MC_inverse(A, P, cdf, n_walks, max_step)

inv_A=zeros(size(A));

for k=1:size(A,1)
      parfor walk=1:n_walks
          previous=k;
          W=1;
          inv_A(k,previous)=inv_A(k,previous)+W;
          i=1;
          while i<=max_step
                aux=rand;
                current=min(find(cdf(previous,:)>aux));
                W=W*A(previous,current)/P(previous,current);                         
                if W==0
                   break;
                end
                i=i+1;
                previous=current;
                inv_A(k,current)=inv_A(k,current)+W;
          end
       end
       inv_A(k,:)=inv_A(k,:)./n_walks;
end
end