function [P, cdf]=prob_forward(A, p)

display('Building of transition matrix');

    if (p == 0)
        P=sparse(zeros(size(A)));
        cdf=sparse(zeros(size(A)));
        for i=1:size(A,1)
                Prow_ind=find(A(i,:));
                num_elem=length(find(A(i,:)));
                for j=1:num_elem
                    P(i,Prow_ind(j))=1/num_elem;
                    cdf(i,Prow_ind(j))=j/num_elem;
                end         
        end

    else

        P=sparse(zeros(size(A)));
        cdf=sparse(zeros(size(A)));
        for i=1:size(A,1)
                Prow_ind=find(A(i,:));
                num_elem=length(find(A(i,:)));
                sump=sum(abs(A(i,Prow_ind)).^p);
                for j=1:num_elem
                    P(i,Prow_ind(j))=abs(A(i,Prow_ind(j))).^p/sump;
                    cdf(i,Prow_ind(j))=sum(P(i,Prow_ind(1:j)));
                end
        end

    end

    P=sparse(P);
    cdf=sparse(cdf);
end
