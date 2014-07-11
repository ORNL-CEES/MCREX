function square_flowmovie(qmethod,sol,tt,By,Bx,A,xy,x,y,bound,ftime,fign,avi)
%SQUARE_FLOWMOVIE generates flow movie on square  domain
%   square_flowmovie(qmethod,U,time,By,Bx,A,xy,x,y,bound,.1,998,0);
%   input
%          qmethod    mixed method 
%          U        velocity solution vector
%          time         solution time vector  
%          By         velocity  y-derivative matrix    
%          Bx         velocity x-derivative matrix    
%          A          vector diffusion matrix
%          xy         velocity nodal coordinate vector   
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          bound      boundary vertex vector
%          ftime      controls speed of animation
%          fign       figure number
%          avi        0/1 switch to generate .avi file
%
% pressure solution is assumed to be essentially zero at outflow 
% so streamfunction satisfies zero Neumann condition there
% calls function xxstreambc.m to set boundary values
%   IFISS function: DJS; 27 April 2012.
% Copyright (c) 2009 D.J. Silvester, H.C. Elman, A. Ramage 
fprintf('\nGenerating flow field movie ... ')
fprintf('\ndo not adjust the figure ...')

L=max(x); nstep=length(tt);
nvtx=length(xy); nu=2*nvtx; 
Asv=A(1:nvtx,1:nvtx); fzero=zeros(nvtx,1);
%
%% apply boundary conditions to the diffusion matrix
[Abc,fzero]=streambc(Asv,fzero,xy,bound);
%[LA,UA]= lu(Abc); 
fprintf('\n   step      min_phi    max_phi\n')
for k=1:nstep
u=sol(:,k); ttk=tt(k);
f=[By,-Bx]*u;
[fsv]=xxstreambc(Asv,f,xy,bound,ttk);
phi=Abc\fsv; %%UA\(LA\fsv);
fprintf('  %4i   %12.5f   %9.3e\n', k, min(phi),max(phi));
%
%% plot velocity
figure(fign)
set(gcf,'Position',[5,5,600,600]);
[X,Y]=meshgrid(x,y);
ax = [min(x)-.1 max(x)+.1 min(y)-.1 max(y)+.1];
xysol = griddata(xy(:,1),xy(:,2),phi,X,Y);
maxphi=max(max(xysol)); minphi=min(min(xysol));
vneg=[minphi:-minphi/24:0];
vpos=[maxphi/6:maxphi/6:maxphi];
vpospos=[0: maxphi/48:maxphi/12];
   contour(X,Y,xysol,[vneg,vpos,vpospos])
 title(['Streamlines : ',num2str(tt(k),'%8.3f'),' seconds']), 
squarex; axis([-1.1,1.1,-1.1,1.1]),  pause(ftime)
%% save frame
if avi==1, MM(k) =  getframe(gcf,[75 75 570 570]);end
end
fprintf('\nThe end.\n')
fprintf('step %g : time is %9.3e\n',k,tt(k))
%%
gohome
cd plotfiles
if avi==1,movie2avi(MM,'cavityflow.avi')
fprintf('\nMovie created: /plotfiles/cavityflow.avi\n'), end
return
