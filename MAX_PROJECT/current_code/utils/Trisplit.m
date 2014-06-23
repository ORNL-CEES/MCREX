function [M, N]=Trisplit(D,Per)

    Num=size(D,1);
    S=Per'*D*Per;

    hold off
    figure()
    spy(S)
    title('Pattern for laplace matrix after redblack reordering')

    M=triu(S,1);
    M=M+diag(diag(S));
    N=-tril(S,-1);

end