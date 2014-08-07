function [sol, rel_res, var, VAR, DX, NWALKS, tally, count, reject]=SEQ_adjoint(fp, dist, P, cdf, numer, stat)

reject=0;

    if (stat.nchecks == 1)
        [sol, rel_res, var, VAR, DX, NWALKS, tally, count]=SEQ_adjoint1(fp, dist, P, cdf, numer, stat); 

    elseif (stat.nchecks == 2) 
        [sol, rel_res, var, VAR, DX, NWALKS, tally, count, reject]=SEQ_adjoint2(fp, dist, P, cdf, numer, stat);

    end
    
end