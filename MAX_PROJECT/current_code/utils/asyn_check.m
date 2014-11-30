function [answer]=asyn_check(A)

B=abs(A);
if max(eigs(B))<1
    answer=true;
else
    answer=false;
end
