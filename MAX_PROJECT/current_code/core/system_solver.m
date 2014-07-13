function [sol, rel_res, var, VAR, DX, NWALKS, iterations, time]=system_solver(scheme, method, fp, dist, P, cdf, numer, stat)

    if strcmp(scheme, 'SEQ') && strcmp(method, 'forward')
        start=cputime;
        [sol, rel_res, var, VAR, DX, NWALKS, iterations]=SEQ_forward(fp, ...
            P, cdf, numer, stat);
        finish=cputime;

    elseif strcmp(scheme, 'MCSA') && strcmp(method, 'forward')
        start=cputime;
        [sol, rel_res, var, VAR, DX, NWALKS, iterations]=MCSA_forward(fp, ...
            P, cdf, numer, stat);
        finish=cputime;
        
    elseif strcmp(scheme, 'SEQ') && strcmp(method, 'adjoint')
        start=cputime;
        [sol, rel_res, var, VAR, DX, NWALKS, iterations]=SEQ_adjoint(fp, ...
            dist, P, cdf, numer, stat);
        finish=cputime;
    

    elseif strcmp(scheme, 'MCSA') && strcmp(method, 'adjoint')
        start=cputime;
        [sol, rel_res, var, VAR, DX, NWALKS, iterations]=MCSA_adjoint(fp, ...
            dist, P, cdf, numer, stat);
        finish=cputime;
           
    else
        error(strcat('Invalid Monte Carlo scheme inserted', '\n', '\n'));
        fprintf('\n\n');
    end
        
    if strcmp(fp.precond, 'gs') || strcmp(fp.precond, 'trisplit') || strcmp(fp.precond, 'alternating')
        sol=(fp.Per)'\sol;
    end
    
    time.start=start;
    time.finish=finish;
end