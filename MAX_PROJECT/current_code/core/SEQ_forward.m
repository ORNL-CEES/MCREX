function [sol, rel_res, VAR, DX, NWALKS, tally, count, reject]=SEQ_forward(fp, P, cdf, numer, stat)

reject=0;

    if (stat.nchecks == 1)
        [sol, rel_res, VAR, DX, NWALKS, tally, count]=SEQ_forward1(fp, P, cdf, numer, stat); 

    elseif (stat.nchecks == 2) 
        [sol, rel_res, VAR, DX, NWALKS, tally, count, reject]=SEQ_forward2(fp, P, cdf, numer, stat);

    end


end