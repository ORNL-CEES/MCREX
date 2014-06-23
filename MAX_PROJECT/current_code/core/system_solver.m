function [sol, rel_err, var, VAR, iterations, time]=system_solver(scheme, method, fp, dist, P, cdf, rich_it, n_walks, max_step, eps)

    if strcmp(scheme, 'SEQ') && strcmp(method, 'forward')
        start=cputime;
        [sol, rel_err, var, VAR, iterations]=SEQ_forward(fp, ...
            P, cdf, rich_it, n_walks, max_step, eps);
        finish=cputime;

    elseif strcmp(scheme, 'MCSA') && strcmp(method, 'forward')
        start=cputime;
        [sol, rel_err, var, VAR, iterations]=MCSA_forward(fp, ...
            P, cdf, rich_it, n_walks, max_step, eps);
        finish=cputime;
        
    elseif strcmp(scheme, 'SEQ') && strcmp(method, 'adjoint')
        start=cputime;
        [sol, rel_err, var, VAR, iterations]=SEQ_adjoint(fp, ...
            dist, P, cdf, rich_it, n_walks, max_step, eps);
        finish=cputime;
    

    elseif strcmp(scheme, 'MCSA') && strcmp(method, 'adjoint')
        start=cputime;
        [sol, rel_err, var, VAR, iterations]=MCSA_adjoint(fp, ...
            dist, P, cdf, rich_it, n_walks, max_step, eps);
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