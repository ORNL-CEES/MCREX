function flag = diag_dominant(A)

    % check for the strong diagonally dominance
    absDiag = abs(diag(A));
    flag_col = all(sum(abs(A), 2) - absDiag < absDiag); 
    flag_row = all(sum(abs(A), 1)' - absDiag < absDiag);

    flag = flag_row && flag_col;

    %check for the simple diagonally dominance
    if ~flag
        flag_col = all(sum(abs(A), 2) - absDiag <= absDiag); 
        flag_row = all(sum(abs(A), 1)' - absDiag <= absDiag);  
        aux = flag_row && flag_col;
        if aux
            aux_row=false;
            aux_col=false;
            for i=1:size(A,1)
                sum_row=sum(abs(A), 1)';
                sum_col=sum(abs(A), 2);
                if sum_row(i) - absDiag(i) <= absDiag(i)
                    aux_row=true;
                end
                if sum_col(i) - absDiag(i) <= absDiag(i)
                    aux_col=true;
                end            
            end
            flag = aux_row && aux_col;
        end

    end
    
end