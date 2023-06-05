function [u, t] = lbm_seqPJ(nx, ny, it)
% lattice Boltzmann example, scalar
% slower, but more didactical version
%  nx,ny  number of lattice points in x,y direction (eg. 33, 257)
%  it     number of iterations (eg. 10000, 350000)
tic;
if nargin < 3
  nx = 33;          % lattice size
  ny = nx;
  it = 10000;       % no. iterations  
end

ux0 = 0.1;          % driving velocity (only in x direction)
Re = 1000;          % Reynolds number

% auxiliary values
q = 9;                   % number of discrete velocities
cs = 1/sqrt(3);
cs2 = 1/cs^2;
cs4 = 1/cs^4;
nu = (ny - 1)*ux0/Re;    % viscosity
tau = 3*nu + 0.5;      % relaxation time of BGK collision operator
om = 1/tau;
omp = 1 - om;

% velocities: D2Q9
% order of velocities (cf. [Krüger et al., table 3.3])
%    7 3 6
%    4 1 2
%    8 5 9
w = [4/9; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36];
cx = [0; 1; 0; -1; 0; 1; -1; -1; 1];
cy = [0; 0; 1; 0; -1; 1; 1; -1; -1];

% careful: x to the right -> is column index, ie 2nd index!
f = ones(ny, nx, q);
feq = zeros(ny, nx, q);

% initialise matrices
% bulk: rho = 1, u = 0 -> feq_i = w_i
f = reshape(f, [nx*ny, q]);
f = f*diag(w);
f = reshape(f, [ny, nx, q]);
% top line: in iteration

for nr=1:it
  % computation of macroscopic variables rho and u
  rho = sum(f, 3);
  ux = (f(:,:,2) + f(:,:,6) + f(:,:,9) ...
        - f(:,:,4) - f(:,:,7) - f(:,:,8))./rho;
  uy = (f(:,:,3) + f(:,:,6) + f(:,:,7) ...
        - f(:,:,5) - f(:,:,8) - f(:,:,9))./rho;
  % top line
  ux(ny,:) = ux0;
  uy(ny,:) = 0;
      
  % computation of the equilibrium distribution feq
  for I=1:q
    feq(:,:,I) = w(I)*rho.*(1 + cs2*(ux*cx(I) + uy*cy(I)) ...
         + 0.5*cs4*(ux*cx(I) + uy*cy(I)).^2 - 0.5*cs2*(ux.*ux + uy.*uy));
  end

  % collision
  % internal points
  f(2:ny-1,2:nx-1,:) = omp*f(2:ny-1,2:nx-1,:) + om*feq(2:ny-1,2:nx-1,:);
  % top line
  f(ny,:,:) = feq(ny,:,:);
  % walls: fullway bounce-back rule
  f(2:ny-1,1,:)  = f(2:ny-1,1,[1,4,5,2,3,8,9,6,7]);   % left
  f(2:ny-1,nx,:) = f(2:ny-1,nx,[1,4,5,2,3,8,9,6,7]);  % right
  f(1,1:nx,:)    = f(1,1:nx,[1,4,5,2,3,8,9,6,7]);     % bottom

  % streaming
  f(:,2:nx,2) = f(:,1:nx-1,2);
  f(2:ny,:,3) = f(1:ny-1,:,3);
  f(:,1:nx-1,4) = f(:,2:nx,4);
  f(1:ny-1,:,5) = f(2:ny,:,5);
  f(2:ny,2:nx,6) = f(1:ny-1,1:nx-1,6);
  f(2:ny,1:nx-1,7) = f(1:ny-1,2:nx,7);
  f(1:ny-1,1:nx-1,8) = f(2:ny,2:nx,8);
  f(1:ny-1,2:nx,9) = f(2:ny,1:nx-1,9);
end

u = sqrt(ux.^2 + uy.^2)/ux0;
% order x/y: C(x,y), origin u(1,1) = left lower corner
u = flipud(u);
t = toc;
