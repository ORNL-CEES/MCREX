function [answer]=MC_converge(H,P)

H_star=zeros(size(H));

for i=1:size(H,1)
    for j=1:size(H,2)
        if (P(i,j)>0)
            H_star(i,j)=sparse((H(i,j).^2)./P(i,j));
        end
    end
end
     
H_star=sparse(H_star);
answer=max(abs(eigs(H_star)))<1;

end