
function [Pc]=prob_adjoint_complement(A, p)
   display('Building of transition matrix');

 
   if (p == 0)    
        P=zeros(size(A));
        for i=1:size(P,1)
                if sum(abs(A(:,i)))>0
                    aux=find(A(:,i));
                    for j=1:length(aux)
                        P(i,aux(j))=(A(aux(j),i)~=0)/(length(aux));
                    end
                end
        end

   else
        P=zeros(size(A));
        for i=1:size(P,1)
                aux=find(A(:,i));
                if sum(abs(A(aux,i)))>0
                    for j=1:length(aux)
                        P(i,aux(j))=abs(A(aux(j),i)).^p/(sum(abs(A(aux,i)).^p));
                    end
                end
        end
   end
    
   
    P=sparse(P);
    
    Pc=zeros(size(A));
    for i=1:size(P,1)
        aux=find(P(i,:));
        l=length(aux);
        for j=1:l
            Pc(i,aux(j))=(1./P(i,aux(j)))/sum(1./P(i,aux));
        end
    end
    
end