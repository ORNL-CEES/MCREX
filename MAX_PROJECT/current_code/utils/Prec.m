function [SET, fp]=Prec(D, prec, fp, SET)

    if strcmp(prec, 'diag')
        fp.prec=diag(diag(D));
        
    elseif strcmp(prec, 'gs')
        [SET.M, SET.N]=GSsplitting(D, fp.Per);
        
    elseif strcmp(prec, 'triblock')
        [SET.M, SET.N]=Triblock(D, fp.Per);
        
     elseif strcmp(prec, 'trisplit')
        [SET.M, SET.N]=Trisplit(D, fp.Per);   
        
    end

end