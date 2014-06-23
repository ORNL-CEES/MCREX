function [u, D, rhs]=laplace(R, G, reac)

    D = delsq(G);
    D = D + reac * eye(size(D));
    % Number of interior points
    Num = sum(G(:)>0);


    %Dirichlet boundary conditions
    rhs = ones(Num,1);
    %rhs=sort(rhs, 'descend');
    if (R == 'N') % For nested dissection, turn off minimum degree ordering.
        spparms('autommd',0)
        u = D\rhs;
        spparms('autommd',1)
    else
        u = D\rhs; % This is used for R=='L' as in this example
    end
    
end


