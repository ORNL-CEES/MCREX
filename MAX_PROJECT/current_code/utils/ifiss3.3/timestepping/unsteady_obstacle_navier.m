%UNSTEADY_OBSTACLE_NAVIER unsteady flow problem in obstacle domain
%   IFISS scriptfile: DJS; 4 November 2013.
% Copyright (c) 2009 D.J. Silvester
clear variables
global pde 
fprintf('Unsteady flow in an obstacle domain ...\n')
global viscosity
viscosity=default('viscosity parameter (default 1/50)',1/50);
%
%% initialize for time stepping iteration: compute Stokes solution
%% load assembled matrices
gohome
cd datafiles
load obstacle_stokes_nobc.mat
%
%
%% set parameters
pde=14; domain=4;
fprintf('Discrete Saddle-Point DAE system ...\n')
tfinal = default('target time? (default 100)',100);
tol = default('accuracy tolerance? (default 1e-4)',1e-4);
nstar = default('averaging frequency? (default 10)',10);
vswitch = default('plot solution evolution? 1/0',0);
dtzero=1e-9; uzero=zeros(size(f));gzero=zeros(size(g));
%% compute solution
tstart = tic;
if qmethod>1, 
np=length(g); 
AxB='defaultAxB'; 
[DT,U,Udot,time] = stabtrNS(qmethod,xy,mv,bound,A,B,sparse(np,np),G,AxB,...
							uzero,dtzero,tfinal,tol,nstar,1);
else  
% call the specialized saddle point system solver 
AxB='colamdAxB'; 
% fixed beta (independent of the viscosity)    
beta=1;   
if qmethod==1,   
	  beta=1/4;     % default parameter
end  
% scaled beta (mutiplied by the viscosity)   
beta=beta*viscosity;
[DT,U,Udot,time] = stabtrNS(qmethod,xy,ev,bound,A,B,beta*C,G,AxB,...
							uzero,dtzero,tfinal,tol,nstar,1);   
end
etoc=toc(tstart); fprintf('Integration took  %8.3e seconds\n\n',etoc) 

save obstacle_unsteadyflow.mat DT U Udot time viscosity
%
% visualise solution and hold the timestep sequence plot
marker = 'bo';
dttplot(DT,99,marker)
if vswitch==1
pause(3)
   if qmethod>1,
   obstacle_unsteadyflowplot(qmethod,mv,U,time,By,Bx,G,xy,xyp,x,y,...
                             bound,bndxy,bnde,obs,0.1,167);
   else
   obstacle_unsteadyflowplot(qmethod,ev,U,time,By,Bx,G,xy,xyp,x,y,...
                             bound,bndxy,bnde,obs,0.1,167);
   end
   fprintf('... to generate a movie run obstacle_flowmovie\n')
   fprintf('    see help obstacle_flowmovie ...\n')
   help obstacle_flowmovie
end
%%
%%% compute divergence error
if qmethod>1
   error_div = q2div(xy,mv,[U(:,end);gzero]);	
else
   error_div = q1div(xy,ev,[U(:,end);gzero]);	 	
end
return
