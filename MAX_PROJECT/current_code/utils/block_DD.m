function [answer, D, H]=block_DD(A, size_block, pol)

    if (mod(size(A,1),size_block)~=0)
            err('the size of the block is not a divisor of the size of the matrix');
    end
    
    n_blocks=size(A,1)/size_block;
    D=zeros(size(A));
    
    for i=1:n_blocks      
        D((i-1)*size_block+1:i*size_block,(i-1)*size_block+1:i*size_block)=A((i-1)*size_block+1:i*size_block,(i-1)*size_block+1:i*size_block);
    end
        
    D=sparse(D);
    
    if strcmp(pol, 'inf')
        H=speye(size(A))-D\A;
        H_tilde=zeros(n_blocks,n_blocks);
        
        for i=1:n_blocks
            for j=1:n_blocks
                if (i~=j)
                    H_tilde(i,j)=norm(H((i-1)*size_block+1:i*size_block,(j-1)*size_block+1:j*size_block),inf);
                end
            end
        end

        answer=norm(H_tilde,inf);
        
    elseif strcmp(pol, '1')
        H=speye(size(A))-A*inv(D);
        H_tilde=zeros(n_blocks,n_blocks);
        for i=1:n_blocks
            for j=1:n_blocks
                if (i~=j)
                    H_tilde(i,j)=norm(H((i-1)*size_block+1:i*size_block,(j-1)*size_block+1:j*size_block),1);
                end
            end
        end
        
        answer=norm(H_tilde,1);
        
     elseif strcmp(pol, '2')
        H=speye(size(A))-D\A;
        H_tilde=zeros(n_blocks,n_blocks); 
        for i=1:n_blocks
            for j=1:n_blocks
                if (i~=j)
                    H_tilde(i,j)=svds(H((i-1)*size_block+1:i*size_block,(j-1)*size_block+1:j*size_block),1);
                end
            end
        end
        
        answer=norm(H_tilde,2);
        
    end
    
return