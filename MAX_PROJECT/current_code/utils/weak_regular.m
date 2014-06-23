function [answer]=weak_regular(B,C)

answer=false;
check1=all( inv(B)>=0 );
check2=all( B\C>=0 );
if min(check1)==1 && min(check2)==1
    answer=true;
end

end