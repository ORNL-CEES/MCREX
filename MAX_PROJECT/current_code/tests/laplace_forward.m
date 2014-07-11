addpath('../core')
addpath('../utils')

%parallel
%poolobj=parpool('local');

shape = 'S'; % Other possible shapes include L,S,N,C,D,A,H,B
n=32;

% creation of te grid
G=gridgen(shape, n);

reac=0;
% computation of the solution with \
[u, D, rhs]=laplace(shape, G, reac);

%algebraic splitting: 'diag', 'gs', 'triblock', 'trisplit', 'alternating'
precond='diag';

%Sequential Monte Carlo or Monte Carlo Synthetic Acceleration
%possible choices: 'SEQ', 'MCSA'
scheme='MCSA';
method='forward';

[fixed_point]=iteration_matrix(precond, D, u, rhs, G);

spy_matrices(fixed_point);

%% Numerical setting

numer.eps=10^(-3);
numer.rich_it=300;%maximal number of Richardson iterations

%% Statistical setting

stat.nwalks=2;
stat.max_step=20;
stat.varcut=0.5;
stat.adapt=1;
dist=1;

%% Definition of the transitional probability

if ~ strcmp(precond, 'alternating')
    [P, cdf]=prob_forward(fixed_point.H, dist);
    
else  
    [P.P1, cdf.cdf1]=prob_forward(fixed_point.H1, dist);
    [P.P2, cdf.cdf2]=prob_forward(fixed_point.H2, dist);
end

%% Solver call

[sol, rel_residual, var, VAR, DX, NWALKS, iterations, time]=system_solver(scheme, method, fixed_point, dist, P, cdf, numer, stat);

%delete(poolobj)

conf=0.025; % error of type 1

plot_sol(G, u, 'ex');
plot_sol(G, sol, scheme);

figure()
plot(sol, '*')
hold on
plot(sol+var*norminv(1-conf/2, 0, 1), 'g*')
plot(sol-var*norminv(1-conf/2, 0, 1), 'g*')
plot(u, 'r*')

count=0;

for i=1:size(u,1)
    if u(i)>sol(i)-var(i)*norminv(1-conf/2, 0, 1) && u(i)<sol(i)+var(i)*norminv(1-conf/2, 0, 1)
        count=count+1;
    end
end

amplitude=0;
for i=1:size(u,1)
    amplitude=amplitude+2*var(i)*norminv(1-conf/2, 0, 1);
end
amplitude=amplitude/norm(u,2);

hold off
bar(NWALKS)


if reac==0
    if strcmp(scheme, 'SEQ')
        save(strcat('../results/diffusion_problem/Sequential_MC/laplace_forward_', scheme, '_', dist, '_', precond));
    elseif strcmp(scheme, 'MCSA')
        save(strcat('../results/diffusion_problem/MCSA/laplace_forward_', scheme, '_', dist, '_', precond));
    end
else
    if strcmp(scheme, 'SEQ')
        save(strcat('../results/DR_problem/Sequential_MC/laplace_forward_', scheme, '_', dist, '_', precond));
    elseif strcmp(scheme, 'MCSA')
        save(strcat('../results/DR_problem/MCSA/laplace_forward_', scheme, '_', dist, '_', precond));
    end      
end
