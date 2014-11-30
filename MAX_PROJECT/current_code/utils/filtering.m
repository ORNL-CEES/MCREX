function [A_tilde]=filtering(A, tol)

    A_tilde=A;
    %norma=norm(A_tilde,'fro');
    maximum=max(abs(A_tilde));
    for i=1:size(A_tilde,1)
        index= abs(A_tilde(i,:)) < tol*maximum;
        A_tilde(i,index)=0;
    end
    
end
