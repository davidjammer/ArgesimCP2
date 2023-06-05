function [uAll, p, tE] = pde_dist(dt, N)
  % partial differential equation example
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
  % split y and delete fixed end points
  u0 = y0(2:p.N); 
  ut0 = y0(N+3:end-1); 
  
  % prepare matrix for ode
  nuk = (p.nu*p.nu)/(p.k*p.k);
  A = diag(ones(N-2,1), -1) - 2*diag(ones(N-1, 1)) + diag(ones(N-2,1), 1);
  A = nuk*A;
  
  spmd
    A = codistributed(A);
    t = 0;
    u = codistributed(u0, codistributor1d(1));
    ut = codistributed(ut0, codistributor1d(1));
    uAll = zeros(p.nSteps + 1, N+1, "codistributed");  % collects u values
    uAll(1,2:N) = u0(1:N-1,:)';
    
    for step=1:p.nSteps
      [t, u, ut] = rk4_step1(@pdeOde, t, dt, u, ut, A);
      uAll(step+1,2:N) = u(1:N-1,:)';
    end
  end
  uAll = gather(uAll);
  tE = toc;

  spmd, cu = getCodistributor(u); cut = getCodistributor(ut); end
  cu = cu{1}; cut = cut{1};
  fprintf("\nPartitioning of u: ");
  fprintf(" %d ", cu.Partition);
  fprintf("\nPartitioning of ut: ");
  fprintf(" %d ", cut.Partition);
  fprintf("\n");
end

function [dx, dv] = pdeOde(t, x, v, A)
  % DGL of the discretized PDE
  dx = v;
  dv = A*x;
end

function [t, xNew, vNew] = rk4_step1(f, t0, h, x, v, A)
  t = t0;
  [k1x, k1v] = f(t, x, v, A);
  [k2x, k2v] = f(t + h/2, x + h*k1x/2, v + h*k1v/2, A);
  [k3x, k3v] = f(t + h/2, x + h*k2x/2, v + h*k2v/2, A);
  [k4x, k4v] = f(t + h/2, x + h*k3x, v + h*k3v, A);

  xNew = x + h*(k1x + 2*k2x + 2*k3x + k4x)/6;
  vNew = v + h*(k1v + 2*k2v + 2*k3v + k4v)/6;
  t = t + h;
end
