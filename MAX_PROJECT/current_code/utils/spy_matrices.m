function []=spy_matrices(fp)

if isfield(fp, 'D')
    figure()
    spy(fp.D)
    title('The 5-point Laplacian')
end

if isfield(fp, 'H')
    figure()
    spy(fp.H)
    title('Iteration matrix H')
end

if isfield(fp, 'M')
    figure()
    spy(fp.M)
    title('Matrix M')
end

if isfield(fp, 'N')
    figure()
    spy(fp.N)
    title('Matrix N')
end 

end