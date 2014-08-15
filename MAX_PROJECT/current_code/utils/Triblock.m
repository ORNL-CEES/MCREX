function [M, N] = Triblock(D, Per)

    Num=size(D,1);
    S=Per'*D*Per;

    hold off
    figure()
    spy(S)
    title('Pattern for laplace matrix after redblack reordering')


    A=S(1:(Num/2), 1:(Num/2));
    B=S(1:(Num/2),(Num/2+1):Num);
    C=S((Num/2+1:Num), 1:(Num/2));
    E=S((Num/2+1:Num),(Num/2+1:Num));

    Schur=E-C*(A\B);
    M=[A zeros(size(A)); C Schur];
    N=[eye(size(A,1)) A\B; zeros(size(A)), eye(size(A,1))];

end