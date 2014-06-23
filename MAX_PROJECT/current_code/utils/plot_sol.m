function []=plot_sol(G, u, scheme)

% Map the solution onto the grid and show it as a contour map.

    figure()
    U = G;
    U(G>0) = full(u(G(G>0)));
    clabel(contour(U));
    prism
    axis square ij
    hold off

    % Now show the solution as a mesh plot.

    n=size(G,1);
    colormap((cool+1)/2);
    mesh(U)
    axis([0 n 0 n 0 max(max(U))])
    if strcmp(scheme, 'ex')
        title('correct solution');
    elseif strcmp(scheme, 'SEQ')
        title('SEQ Monte Carlo solution');
    elseif strcmp(scheme, 'MCSA')
        title('MCSA Monte Carlo solution');
    end
    axis square ij

end