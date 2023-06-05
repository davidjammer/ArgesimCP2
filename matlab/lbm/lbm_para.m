function [u,t] = lbm_para(Nx, Ny, REP)
tic;

RHO_0 = 1.0;
u_0 = 0.1;
RE = 1000;


%0 = Fluid Cell
%1 = Wall Cell
%2 = Driving Cell

geometry = ones(Ny,Nx);
geometry(2:end-1,2:end-1) = 0;
geometry(1,:) = 2;

C = 1; 
E = 2;  S = 3; 
W = 4;  N = 5;
NE = 6; SE = 7;
SW = 8; NW = 9;


spmd
	%calculate the viscosity [Equation 6]
	viscosity = (Ny - 1) * u_0 / RE;
	%calculate tau [Equation 5]
	tau = (6 * viscosity + 1) / 2;
    
	
	numx = 1;
	numy = 1;

	while (numx * numy) < numlabs
		if numx == numy
			numx = numx * 2;
		else
			numy = numy * 2;
		end
	end

	nx = Nx / numx;
	ny = Ny / numy;
	
	rank = [fix((labindex-1)/numx + 1), fix(mod((labindex-1), numx) + 1)];
	RANK = reshape(1:numlabs, numx, numy)';
	index = [(rank(1) - 1) * ny, (rank(2) - 1) * nx];
	
	l_geometry = geometry(index(1)+1:index(1)+ny, index(2)+1:index(2)+nx);
	
	[driving_IDs.y,driving_IDs.x] = find(l_geometry == 2);
	[wall_IDs.y,wall_IDs.x] = find(l_geometry == 1);
	[fluid_IDs.y,fluid_IDs.x] = find(l_geometry == 0);
	 
	
	rho = zeros(ny,nx);
	ux = zeros(ny,nx);
	uy = zeros(ny,nx);
	usqr = zeros(ny,nx);
	
	f = zeros(ny,nx,9);
	feq = zeros(ny,nx,9);

	%init dist func values
	f(:,:,C) = RHO_0 *4.0 / 9.0;
	f(:,:,N) = RHO_0 / 9.0;
	f(:,:,E) = RHO_0 / 9.0;
	f(:,:,S) = RHO_0 / 9.0;
	f(:,:,W) = RHO_0 / 9.0;
	f(:,:,NE) = RHO_0 / 36.0;
	f(:,:,SE) = RHO_0 / 36.0;
	f(:,:,SW) = RHO_0 / 36.0;
	f(:,:,NW) = RHO_0 / 36.0;
	
	%sim loop
	for rep=1:REP
        disp(rep);
	   %collision step
	  
	   rho = sum(f,3);
	   ux = (f(:,:,E) - f(:,:,W) + f(:,:,NE) + f(:,:,SE) - f(:,:,SW) - f(:,:,NW)) ./rho;
	   uy = (f(:,:,N) - f(:,:,S) + f(:,:,NE) + f(:,:,NW) - f(:,:,SE) - f(:,:,SW)) ./rho; 
	   
	   ux(driving_IDs.y,driving_IDs.x) = u_0;
	   uy(driving_IDs.y,driving_IDs.x) = 0;
	   usqr(:) = ux.*ux + uy.*uy;

	   
	   
	   feq(:,:,C) = 4.0 / 9.0 * rho .* (1 - 1.5 * usqr);
	   
	   feq(:,:,N) = 1.0 / 9.0 * rho .* (1.0 + 3.0 * uy + 4.5 * uy.^2 - 1.5 * usqr);
	   feq(:,:,E) = 1.0 / 9.0 * rho .* (1.0 + 3.0 * ux + 4.5 * ux.^2 - 1.5 * usqr);
	   feq(:,:,S) = 1.0 / 9.0 * rho .* (1.0 - 3.0 * uy + 4.5 * uy.^2 - 1.5 * usqr);
	   feq(:,:,W) = 1.0 / 9.0 * rho .* (1.0 - 3.0 * ux + 4.5 * ux.^2 - 1.5 * usqr);
	   
	   feq(:,:,NE) = 1.0 / 36.0 * rho .* (1.0 + 3.0 * ( ux + uy) + 4.5 * ( ux + uy).^2 - 1.5 * usqr);
	   feq(:,:,SE) = 1.0 / 36.0 * rho .* (1.0 + 3.0 * ( ux - uy) + 4.5 * ( ux - uy).^2 - 1.5 * usqr);
	   feq(:,:,SW) = 1.0 / 36.0 * rho .* (1.0 + 3.0 * (-ux - uy) + 4.5 * (-ux - uy).^2 - 1.5 * usqr);
	   feq(:,:,NW) = 1.0 / 36.0 * rho .* (1.0 + 3.0 * (-ux + uy) + 4.5 * (-ux + uy).^2 - 1.5 * usqr);

	   
	  
	   for x=1:nx
		  for y=1:ny
			
			 if l_geometry(y,x) == 0 %%fluid
				f(y,x,:) = f(y,x,:) * (1.0 - 1.0 / tau) + feq(y,x,:) /tau;			
			 elseif l_geometry(y,x) == 2 %%driving
				f(y,x,:) = feq(y,x,:);			
			 else %%wall
				f(y,x,[E W]) = f(y,x, [W E]);
				f(y,x,[N S]) = f(y,x, [S N]);
				f(y,x,[NE SW]) = f(y,x, [SW NE]);
				f(y,x,[SE NW]) = f(y,x, [NW SE]);
			 end
		  end
	   end
	   

	   %%Osten
	   labBarrier;
	   f(:,:,E) = shift_right(f(:,:,E),rank,RANK);
	   %%Norden
	   labBarrier;
	   f(:,:,N) = shift_up(f(:,:,N),rank,RANK);
	   %%Sueden
	   labBarrier;
	   f(:,:,S) = shift_down(f(:,:,S),rank,RANK);
	   %%Westen
	   labBarrier;
	   f(:,:,W) = shift_left(f(:,:,W),rank,RANK);
	   %%NE
 	   labBarrier;
 	   f(:,:,NE) = shift_right(f(:,:,NE),rank,RANK);
 	   labBarrier;
 	   f(:,:,NE) = shift_up(f(:,:,NE),rank,RANK);
	   %%SE
	   labBarrier;
	   f(:,:,SE) = shift_right(f(:,:,SE),rank,RANK);
	   labBarrier;
	   f(:,:,SE) = shift_down(f(:,:,SE),rank,RANK);
	   %%SW
	   labBarrier;
	   f(:,:,SW) = shift_left(f(:,:,SW),rank,RANK);
	   labBarrier;
	   f(:,:,SW) = shift_down(f(:,:,SW),rank,RANK);
	   %%NW
	   labBarrier;
	   f(:,:,NW) = shift_left(f(:,:,NW),rank,RANK);
	   labBarrier;
	   f(:,:,NW) = shift_up(f(:,:,NW),rank,RANK);
	end
	
end

n = length(ux);
UX = zeros(Ny,Nx);
UY = zeros(Ny,Nx);

for i=1:n
    INDEX = index{i};
    uxl = ux{i};
    uyl = uy{i};
    UX(INDEX(1)+1:INDEX(1)+ny{i}, INDEX(2)+1:INDEX(2)+nx{i}) = uxl(1:ny{i},1:nx{i});
    UY(INDEX(1)+1:INDEX(1)+ny{i}, INDEX(2)+1:INDEX(2)+nx{i}) = uyl(1:ny{i},1:nx{i});
end

u=sqrt(UX.^2+UY.^2)/u_0; % calculate relative macroscopic velocity magnitude


t = toc;


end