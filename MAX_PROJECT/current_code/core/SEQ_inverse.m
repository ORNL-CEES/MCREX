function [inv_A]=SEQ_inverse(A, n_walks, max_step, dist)

diagonally_dominant=1;

for i=1:size(A,1)
    if abs(A(i,i))<(sum(abs(A(i,:)))-abs(A(i,i)))
        diagonally_dominant=0;
    end
end

    
if (diagonally_dominant==1) && (abs(eig(A))<1)
    M=eye(size(A))-A;
    [P,cdf]=prob_forward(M, dist);
    [inv_A]=MC_inverse(M, P, cdf, n_walks, max_step);
end

if (diagonally_dominant==1) && (abs(eig(A))>=1)
    D=diag(diag(A));
    inv_D=diag(1./diag(D));
    T=eye(size(A))-D\A;
    Y=eye(size(A))-T;
    M=eye(size(Y))-Y\inv_D;
    [P,cdf]=prob_forward(M, dist);
    [inv_A]=MC_inverse(M, P, cdf, n_walks, max_step);
end

if (diagonally_dominant==0)
    err('code for generic matrix still needs to be implemented');
end

end