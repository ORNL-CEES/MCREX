dimen=50;
A=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
rhs=[1:50]';
u=A\rhs;


% a priori knowledge
A1=2*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
A2=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),-1);
A3=4*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1);
A4=4*diag(ones(dimen,1)) - diag(rand(dimen-1,1),1) - diag(rand(dimen-1,1),-1);

Prec=diag(diag(A));
Prec1=diag(diag(A1));
Prec4=diag(diag(A4));

H=eye(size(A))-1/1*inv(Prec)*A;

H1=eye(size(A1))-1/1*inv(Prec1)*A1;
H2=eye(size(A2))-1/2*A2;
H3=eye(size(A3))-1/2*A3;
H4=eye(size(A4))-1/1*inv(Prec4)*A4;


rhs=1/1*inv(Prec)*rhs;

P=zeros(size(H));
for i=1:size(P,1)
    for j=1:size(P,2)
        P(i,j)=abs(H(i,j))/(sum(abs(H(i,:))));
    end
end

P1=zeros(size(H));
for i=1:size(P1,1)
    for j=1:size(P1,2)
        P1(i,j)=abs(H1(i,j))/(sum(abs(H1(i,:))));
    end
end

P2=zeros(size(H));
for i=1:size(P2,1)
    for j=1:size(P2,2)
        P2(i,j)=abs(H2(i,j))/(sum(abs(H2(i,:))));
    end
end

P3=zeros(size(H));
for i=1:size(P3,1)
    for j=1:size(P3,2)
        P3(i,j)=abs(H3(i,j))/(sum(abs(H3(i,:))));
    end
end

P4=zeros(size(H));
for i=1:size(P4,1)
    for j=1:size(P4,2)
        P4(i,j)=abs(H4(i,j))/(sum(abs(H4(i,:))));
    end
end

%update the distribution P
Pbay=zeros(size(H));
for i=1:size(Pbay,1)
    for j=1:size(Pbay,2)
        if (P1(i,j)+P2(i,j)+P3(i,j)+P4(i,j))~=0
            Pbay(i,j)=(P(i,j)+P1(i,j)+P2(i,j)+P3(i,j)+P4(i,j)/5)*(P(i,j))/((P1(i,j)+P2(i,j)+P3(i,j)+P4(i,j))/4);
        end
    end
    Pbay(i,:)=Pbay(i,:)./norm(Pbay(i,:),1);
end

eps=10^(-3);
n_walks=10000;

cdf=Pbay;
%computation of the cumulative probability
for i=1:size(cdf,1)
    for j=2:size(cdf,2)
        aux=max(find(cdf(i,1:j-1)))+0;
        if (cdf(i,j)~=0 && ~isempty(aux))
            cdf(i,j)=cdf(i,j)+cdf(i,aux);
        end
    end
end

max_step=10;

%%


u_approx=MC_forward(H, rhs, Pbay, cdf, n_walks, max_step);
rel_error=sqrt(sum((u-u_approx).^2))/sqrt(sum((u).^2));

