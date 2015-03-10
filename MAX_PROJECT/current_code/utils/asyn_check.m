function [answer]=asyn_check(A,b)

B=abs(A);
if max(eigs(B))<1
    answer=MC_converge(A,b);
end
