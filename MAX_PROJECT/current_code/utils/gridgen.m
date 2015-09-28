function [G]=gridgen(shape, n)

    % Generate and display the grid.
    G = numgrid(shape,n);
    %figure()
    %spy(G)
    %title('A finite difference grid');
    % Show a smaller version as sample.
    %g = numgrid(R,12)

end

