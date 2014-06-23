function [answer]=asyn_check(A)

B=abs(A);
if max(eig(B))<1
    answer=true;
else
    answer=false;
end
