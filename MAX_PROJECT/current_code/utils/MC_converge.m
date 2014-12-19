function [answer]=MC_converge(H,P)

display('MC_check');

H_star=zeros(size(H));

for i=1:size(H,1)
    aux=find(H(:,i));
    for j=1:length(aux)
        if (P(i,aux(j))>0)
            H_star(i,aux(j))=sparse((H(i,aux(j)).^2)./P(i,aux(j)));
        end
    end
end
     
H_star=sparse(H_star);
answer=max(abs(eigs(H_star)))<1;

end