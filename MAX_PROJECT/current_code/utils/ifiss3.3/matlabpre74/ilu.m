function [L,U] = ilu(A,setup)
%ILU  incomplete ILU with no fill-in   
%   [L,U] = ilu(A,setup);
%   input
%          A          coefficient matrix 
%          setup      factorization structure (not used)
%   output
%          L,U        sparse factors
%
% calls luinc: but is slow compared to matlab ilu!
%   IFISS function: DJS; 31 January 2010.
% Copyright (c) 2009 D.J. Silvester, H.C. Elman, A. Ramage 

[L,U]=luinc(A,'0');
return
