function [sol, rel_res, VAR, REL_RES, DX, NWALKS, tally, count, reject]=MCSA_forward(fp, P, cdf, numer, stat)

reject=0;

    if (stat.nchecks == 1)
        [sol, rel_res, VAR, REL_RES, DX, NWALKS, tally, count]=MCSA_forward1(fp, P, cdf, numer, stat); 

    elseif (stat.nchecks == 2) 
        [sol, rel_res, VAR, REL_RES, DX, NWALKS, tally, count, reject]=MCSA_forward2(fp, P, cdf, numer, stat);

    end


end