function B=diag_reorder(A)

    A2=abs(A);
    p=zeros(size(A,1),1);

    for i=1:size(A2,1)
        p=find(A2(i,:)==max(A2(i,:)));
    end

    B=A(:,p);
end