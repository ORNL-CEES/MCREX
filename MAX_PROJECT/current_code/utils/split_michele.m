function [Z,W,D_inv] = split_michele(A)

    % in input we have the compressed factors for the approximate inverse of
    % A

    Z=triu(A,-1);
    W=tril(A,1);
    D_inv=diag(diag(A));
    Z=Z+speye(size(A,1));
    W=W+speye(size(A,1));
    D_inv=sparse(D_inv);
end