function [answer]=P_regular(B,C)

answer=false;
if min(eig(B+C))>0
    answer=true;
end

end