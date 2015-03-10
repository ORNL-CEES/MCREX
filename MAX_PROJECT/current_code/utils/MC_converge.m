function [answer]=MC_converge(H,P)

display('MC_check');

H_star=zeros(size(H));

for i=1:size(H,1)
    aux=find(H(:,i));
    l=length(aux);
    for j=1:l
        if (P(i,aux(j))>0)
            H_star(i,aux(j))=(H(i,aux(j)).^2)./P(i,aux(j));
        end
    end
end
     
H_star=sparse(H_star);
answer=max(abs(eigs(H_star)));
answer

end