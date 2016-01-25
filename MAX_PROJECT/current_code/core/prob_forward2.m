function [P]=prob_forward2(A, p)

display('Building of transition matrix');

    if (p == 0)
        P=sparse(zeros(size(A)));
        for i=1:size(A,1)
                Prow_ind=find(A(i,:));
                num_elem=length(find(A(i,:)));
                for j=1:num_elem
                    P(i,Prow_ind(j))=1/num_elem;
                end         
        end

    else

        P=sparse(zeros(size(A)));
        for i=1:size(A,1)
                Prow_ind=find(A(i,:));
                num_elem=length(Prow_ind);
                sump=sum(abs(A(i,Prow_ind)).^p);
                for j=1:num_elem
                    P(i,Prow_ind(j))=abs(A(i,Prow_ind(j))).^p/sump;
                end
        end

    end

    P=sparse(P);
    
end
