function [P]=prob_adjoint3(A, p)

    display('Building of transition matrix');

 
   if (p == 0)    
        P=sparse(zeros(size(A)));
        for i=1:size(A,2)    
            aux=find(A(:,i));
            l=length(aux);
            for j=1:l
                    P(i,aux(j))=1/l;
            end           
        end

   else
        P=sparse(zeros(size(A)));
        for i=1:size(A,2)
            aux=find(A(:,i));
            sump=(sum(abs(A(aux,i)).^p));
            l=length(aux);                   
            for j=1:l
                 P(i,aux(j))=abs(A(aux(j),i)).^p/sump;
            end                
        end
   end
    
   
    P=sparse(P);
 
end