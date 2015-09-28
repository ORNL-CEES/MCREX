function [fp]=iteration_matrix2(precond, D, u, rhs, size_block)

fp.precond=precond;
fp.D=D;
MATRICES=struct;
SET=struct;
SET.size_block = size_block;

    [~, fp]=Prec(D, precond, fp, SET);
    fp.H=sparse(eye(size(D,1))-fp.prec\(D));
    fp.rhs=(fp.prec)\rhs;
    fp.u=u;

end