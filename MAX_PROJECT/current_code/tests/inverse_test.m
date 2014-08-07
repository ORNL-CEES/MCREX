addpath('../core')
addpath('../utils')

%parpool('local')

dimen=10;
A=1./16*(8*diag(ones(dimen,1)) - diag(ones(dimen-1,1),1) - diag(ones(dimen-1,1),-1));

max_step=20;
err=[];
WALKS=[];
n_walks=10^4;
dist=1;

for i=1:6
    n_walks=10^i;
    WALKS=[WALKS n_walks];
    [inv_A]=SEQ_inverse(A, n_walks, max_step, dist);
    err=[err norm(inv(A)-inv_A,'fro')/norm(inv(A),'fro')];
end

loglog(WALKS, err, '-or');
hold on
loglog(WALKS, 1./sqrt(WALKS));


% B=rand(dimen,dimen);
% [inv_B]=MCSA_inverse(B, n_walks, max_step);

save('../results/inverse/inverse_test');