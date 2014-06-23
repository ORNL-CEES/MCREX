dimen=10;
A=2*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1);
rhs=ones(dimen,1);
u=A\rhs;

Prec=diag(diag(A));

H=eye(size(A))-1/2*inv(Prec)*A;
rhs=1/2*inv(Prec)*rhs;

P=zeros(size(H));
for i=1:size(P,1)
    for j=1:size(P,2)
        P(i,j)=abs(H(i,j))/(sum(abs(H(i,:))));
    end
end

eps=10^(-6);
n_walks=500;

PW=zeros(n_walks,size(rhs,1));
X=zeros(n_walks,size(rhs,1));

cdf=P;
%computation of the cumulative probability
for i=1:size(cdf,1)
    for j=2:size(cdf,2)
        aux=max(find(cdf(i,1:j-1)))+0;
        if (cdf(i,j)~=0 && ~isempty(aux))
            cdf(i,j)=cdf(i,j)+cdf(i,aux);
        end
    end
end

max_step=1000;
%%

for k=1:size(H,1)
    for walk=1:n_walks
    previous=k;
    PW(walk,k)=1;
    W=1;
    current=k;
    X(walk,k)=X(walk,k)+W*rhs(current);
    i=1;
        while i<=max_step
            aux=rand;
            current=min(find(cdf(previous,:)>aux));
            PW(walk,k)=PW(walk,k)*P(previous,current);
            W=W*H(previous,current)/P(previous,current);
            X(walk,k)=X(walk,k)+W*rhs(current);
            i=i+1;
            previous=current;
        end
    end
end

u_approx=mean(X,1);
