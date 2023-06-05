function tE = pde_seq(dt, N)
  % partial differential equation example, scalar
  % solving by "brute force", i.e. with full matrix
  %   dt          time step size of output values
  %   N           number of intervals for space discretisation
  if nargin < 2
    dt = 0.01;
    N = 500;
  end
  
  tic;
  p = setParams(dt, N);  
  y0 = setInitialCondition(p);     % incl. fixed end points of u and u_t
  y0 = [y0(2:p.N); y0(N+3:end-1)]; % delete fixed end points
  
  % prepare matrix for ode
  nuk = (p.nu*p.nu)/(p.k*p.k);
  A = diag(ones(N-2,1), -1) - 2*diag(ones(N-1, 1)) + diag(ones(N-2,1), 1);
  A = nuk*A;
  fOde = @(t,y) pdeOde(t, y, A);

  t = 0;
  y = y0;
  uAll = zeros(p.nSteps + 1, N+1);   % collects u values
  uAll(1,2:N) = y0(1:N-1,:)';
  for step=1:p.nSteps
    [t, y] = rk4_step(fOde, t, dt, y);
    uAll(step+1,2:N) = y(1:N-1,:)';
  end
  
  tE = toc;
  plotResults(uAll, p);
end

function dy = pdeOde(t, y, A)
  % DGL of the discretized PDE
  n = numel(y)/2;     % number of unknowns

  u = y(1:n);
  ut = y(n+1:end);
  dx = ut;
  dv = A*u;

  dy = [dx; dv];
end
