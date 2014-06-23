function [M, N]=GSsplitting(D, Per)

    Num=size(D,1);
    S=Per'*D*Per;
    D1=S(1:Num/2, 1:Num/2);
    D2=S((Num/2+1):Num, (Num/2+1):Num);
    B=S((Num/2+1):Num, 1:Num/2);

    M=[D1 zeros(size(D1)); B D2];
    N=[zeros(size(D1)) -B'; zeros(size(D1)) zeros(size(D1))];

end