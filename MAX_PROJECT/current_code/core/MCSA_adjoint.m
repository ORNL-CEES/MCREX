function [sol, rel_res, VAR, RES, DX, NWALKS, tally, count, reject]=MCSA_adjoint(fp, dist, P, cdf, numer, stat)

reject=0;

    if (stat.nchecks == 1)
        [sol, rel_res, VAR, RES, DX, NWALKS, tally, count]=MCSA_adjoint1(fp, dist, P, cdf, numer, stat); 

    elseif (stat.nchecks == 2) 
        [sol, rel_res, VAR, RES, DX, NWALKS, tally, count, reject]=MCSA_adjoint2(fp, dist, P, cdf, numer, stat);

    end
    
end