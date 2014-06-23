X=normrnd(0,1,10000,1);

Mean=zeros(size(X));
err=zeros(size(X));

for i=1:size(X,1)
    Mean(i)=mean(X(1:i));
    err(i)=Mean(i)-0;
end

semilogx([1:10000]', Mean);
hold off
figure()
semilogx([1:10000]', err);