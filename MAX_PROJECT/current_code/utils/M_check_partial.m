function flag = M_check_partial(A)

    % N.B.: this is just a partial check of the M property since the
    % semipositive definiteness of the matrix is not checked
    
    Diag=diag(diag(A));
    N=A-Diag;
    flag=all(all(Diag >= 0)~=0);
    if flag
        flag= flag && all(all(N <= 0)~=0);
    end
    
end