function [u, t] = lbm_seq(nx, ny, iterations)
tic;
%2D square Lattice Boltzmann Method reference implementation
%
%Problem description can be found in:
%Hou et al.: Simulation of Cavity Flow by the Lattice Boltzmann Method.
%Journal of Computational Physics, Vol. 118, Issue 2, 05/1995.
%
%Author: Rene Fink, Wismar University of Applied Sciences, Research Group
%Computational Engineering and Automation
%
%Date: 2006-09-26
% begin experiment parameter settings

rho_0=1; % initial density
u_0=0.1; % driving velocity on the top boundary
Re=1000; % Reynolds number

C=1;E=2;S=3;W=4;N=5;NE=6;SE=7;SW=8;NW=9; % helper variables for directions
% directions are indiced as follows:
%   9 5 6
%   4 1 2
%   8 3 7

geometry=ones(ny,nx); % wall cells (geometry==1)
geometry(2:ny-1,2:nx-1)=0; % fluid cells (geometry==0)
geometry(1,:)=2; % driving cells on the top boundary (geometry==2)

viscosity=(ny-1)*u_0/Re; % kinematic viscosity (0.005 <= viscosity <= 0.2)
tau=(6*viscosity+1)/2; % relaxation time

f=zeros(nx*ny,9);     % distribution function values of each cell
feq=zeros(nx*ny,9);   % equilibrium disribution function value
rho=zeros(nx*ny,1);   % macroscopic density
ux=zeros(nx*ny,1);    % macroscopic velocity in x direction
uy=zeros(nx*ny,1);    % macroscopic velocity in y direction
usqr=zeros(nx*ny,1);  % ux^2 + uy^2

% set initial distribution
f(:,C)=rho_0*4/9;
f(:,[E S W N])=rho_0/9;
f(:,[NE SE SW NW])=rho_0/36;

% create indices of all fluid, wall and driving cells
FL=find(geometry==0);
WALL=find(geometry==1);
DR=find(geometry==2);

for i=1:iterations
    % computation of macroscopic variables rho and u
    rho(:)=sum(f,2);
    ux(:)=(f(:,E)-f(:,W)+f(:,NE)+f(:,SE)-f(:,SW)-f(:,NW))./rho;
    uy(:)=(f(:,N)-f(:,S)+f(:,NE)+f(:,NW)-f(:,SE)-f(:,SW))./rho;
    % set velocities for driving cells
    ux(DR)=u_0;
    uy(DR)=0;
    usqr(:)=ux.*ux+uy.*uy;

    % computation of the equilibrium distribution
    feq(:,C)=(4/9)*rho.*(1-1.5*usqr);
    feq(:,E)=(1/9)*rho.*(1+3*ux+4.5*ux.^2-1.5*usqr);
    feq(:,S)=(1/9)*rho.*(1-3*uy+4.5*uy.^2-1.5*usqr);
    feq(:,W)=(1/9)*rho.*(1-3*ux+4.5*ux.^2-1.5*usqr);
    feq(:,N)=(1/9)*rho.*(1+3*uy+4.5*uy.^2-1.5*usqr);
    feq(:,NE)=(1/36)*rho.*(1+3*(ux+uy)+4.5*(ux+uy).^2-1.5*usqr);
    feq(:,SE)=(1/36)*rho.*(1+3*(ux-uy)+4.5*(ux-uy).^2-1.5*usqr);
    feq(:,SW)=(1/36)*rho.*(1+3*(-ux-uy)+4.5*(-ux-uy).^2-1.5*usqr);
    feq(:,NW)=(1/36)*rho.*(1+3*(-ux+uy)+4.5*(-ux+uy).^2-1.5*usqr);
    
    % collision step
    % walls: fullway bounce-back rule
    f(WALL,[C E S W N NE SE SW NW])=f(WALL,[C W N E S SW NW NE SE]);
    % driving cells
    f(DR,:)=feq(DR,:); % distribution function value = equilibrium value
    % fluid cells
    f(FL,:)=f(FL,:)*(1-1/tau)+feq(FL,:)/tau;

    % propagation step
    f=reshape(f,[ny,nx,9]); % transform f for easy propagation

    % particle propagation in directions E,N,S,W,NE,SE,SW,NW
    f(:,2:nx,E)=f(:,1:nx-1,E);
    f(2:ny,:,S)=f(1:ny-1,:,S);
    f(:,1:nx-1,W)=f(:,2:nx,W);
    f(1:ny-1,:,N)=f(2:ny,:,N);
    f(1:ny-1,2:nx,NE)=f(2:ny,1:nx-1,NE);
    f(2:ny,2:nx,SE)=f(1:ny-1,1:nx-1,SE);
    f(2:ny,1:nx-1,SW)=f(1:ny-1,2:nx,SW);
    f(1:ny-1,1:nx-1,NW)=f(2:ny,2:nx,NW);
    
    % re-transform f for next iteration step
    f=reshape(f,[nx*ny,9]);
end

% calculation of relative macroscopic velocity magnitude
u=sqrt(ux.^2+uy.^2)/u_0;
u=reshape(u,ny,nx);
t = toc;

end
