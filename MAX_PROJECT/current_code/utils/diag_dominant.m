function flag = diag_dominant(A)

    % check for the strong diagonally dominance
    absDiag = abs(diag(A));
    flag_row_s = all(sum(abs(A), 2) - absDiag < absDiag); 
    flag_col_s = all(sum(abs(A), 1)' - absDiag < absDiag);

    if flag_row_s
        display('strictly diagonally domimant by rows');
    end
    
    if flag_col_s
        display('strictly diagonally domimant by columns');
    end

    flag = flag_row_s && flag_col_s; 
    
    %check for the simple diagonally dominance
    if ~flag
        flag_row_d = all(sum(abs(A), 2) - absDiag <= absDiag); 
        flag_col_d = all(sum(abs(A), 1)' - absDiag <= absDiag);  
        
        if flag_row_d && ~flag_row_s
            display('weakly diagonally domimant by rows');
        end
            
        if flag_col_d && ~flag_col_s
            display('weakly diagonally domimant by columns');
        end
        
        flag = flag_row_d || flag_col_d; 
    
        if ~flag
            display('NO diagonal dominance');
        end
        
    end
 
end