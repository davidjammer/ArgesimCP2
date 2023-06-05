function [u,t] = lbm_para1(Nx, Ny, iterations)
tic;

rho_0=1; % initial density
u_0=0.1; % driving velocity on the top boundary
Re=1000; % Reynolds number

C=1;E=2;S=3;W=4;N=5;NE=6;SE=7;SW=8;NW=9; % helper variables for directions
% directions are indiced as follows:
%   9 5 6
%   4 1 2
%   8 3 7

geometry=ones(Ny,Nx); % wall cells (geometry==1)
geometry(2:Ny-1,2:Nx-1)=0; % fluid cells (geometry==0)
geometry(1,:)=2; % driving cells on the top boundary (geometry==2)

spmd
    % compute grid of workers assuming numlabs is power of 2
    numx = pow2(ceil(log2(numlabs)/2));
    numy = numlabs/numx;
	nx = Nx / numx;
	ny = Ny / numy;
	
    % 2d-coordinates of worker in core grid
	rank = [fix((labindex-1)/numx + 1), fix(mod((labindex-1), numx) + 1)];
    % position of all workers in core grid
	RANK = reshape(1:numlabs, numx, numy)';
    % position of local array in global array
	index = [(rank(1) - 1) * ny, (rank(2) - 1) * nx]
	
	l_geometry = geometry(index(1)+1:index(1)+ny, index(2)+1:index(2)+nx);

    viscosity=(Ny-1)*u_0/Re; % kinematic viscosity (0.005 <= viscosity <= 0.2)
    tau=(6*viscosity+1)/2;   % relaxation time

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
    FL=find(l_geometry==0);
    WALL=find(l_geometry==1);
    DR=find(l_geometry==2);

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
	   labBarrier;
	   f(:,:,E) = shift_E(f(:,:,E),rank,RANK);
	   labBarrier;
	   f(:,:,N) = shift_N(f(:,:,N),rank,RANK);
	   labBarrier;
	   f(:,:,S) = shift_S(f(:,:,S),rank,RANK);
	   labBarrier;
	   f(:,:,W) = shift_W(f(:,:,W),rank,RANK);
 	   labBarrier;
 	   f(:,:,NE) = shift_NE(f(:,:,NE),rank,RANK);
	   labBarrier;
	   f(:,:,SE) = shift_SE(f(:,:,SE),rank,RANK);
	   labBarrier;
	   f(:,:,SW) = shift_SW(f(:,:,SW),rank,RANK);
	   labBarrier;
	   f(:,:,NW) = shift_NW(f(:,:,NW),rank,RANK);

	   % re-transform f for next iteration step
	   f=reshape(f,[nx*ny,9]); 
    end

end

% collecting local arrays of ux and uy in global ones
n = length(ux);
UX = zeros(Ny,Nx);
UY = zeros(Ny,Nx);
for i=1:n
    INDEX = index{i};
    uxl = reshape(ux{i},ny{i},nx{i});
    uyl = reshape(uy{i},ny{i},nx{i});
    UX(INDEX(1)+1:INDEX(1)+ny{i}, INDEX(2)+1:INDEX(2)+nx{i}) = uxl(1:ny{i},1:nx{i});
    UY(INDEX(1)+1:INDEX(1)+ny{i}, INDEX(2)+1:INDEX(2)+nx{i}) = uyl(1:ny{i},1:nx{i});
end

% calculation of relative macroscopic velocity magnitude
u=sqrt(UX.^2+UY.^2)/u_0;
t = toc;

end