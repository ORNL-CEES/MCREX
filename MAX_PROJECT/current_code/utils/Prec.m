function [SET, fp]=Prec(D, prec, fp, SET)

    if strcmp(prec, 'diag')
        fp.prec=diag(diag(D));
        
    elseif strcmp(prec, 'block_diag')
        size_block=SET.size_block;
        n_blocks=size(D,1)/size_block;

        D_block=sparse(zeros(size(D)));

        for i=0:n_blocks-1
            D_block(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block) = D(i*size_block+1:(i+1)*size_block,i*size_block+1:(i+1)*size_block);
        end

        fp.prec=sparse(D_block);

        
    elseif strcmp(prec, 'gs')
        [SET.M, SET.N]=GSsplitting(D, fp.Per);
        
    elseif strcmp(prec, 'triblock')
        [SET.M, SET.N]=Triblock(D, fp.Per);
        
    elseif strcmp(prec, 'trisplit')
        [SET.M, SET.N]=Trisplit(D, fp.Per); 
        
    elseif strcmp(prec, 'no')   
        fp.prec = speye(size(D));
    end
    
end