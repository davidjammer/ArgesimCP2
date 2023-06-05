function [uAll, p, tE] = pde_distslow(dt, N)
  % partial differential equation example
  % Attention: VERY slow, since y is awkwardly distributed
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
  
  spmd
    A = codistributed(A);
    t = 0;
    y = codistributed(y0, codistributor1d(1)); 
       % wrong distribution, e.g. for n=2: 
       % upper half (=u) at 1, lower half (= ut) at 2
    uAll = zeros(p.nSteps + 1, N+1, "codistributed");  % collects u values
    uAll(1,2:N) = y0(1:N-1,:)';
  
    for step=1:p.nSteps
      [t, y] = rk4_step1(@pdeOde, t, dt, y, A);
      uAll(step+1,2:N) = y(1:N-1,:)';
    end
  end
  uAll = gather(uAll);
  tE = toc;

  spmd, cy = getCodistributor(y); end
  cy = cy{1};
  fprintf("\nPartitioning of y: ");
  fprintf(" %d ", cy.Partition);
  fprintf("\n");
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

function [t, yNew] = rk4_step1(f, t0, h, y, A)
  f1 = @(t,y) f(t, y, A);
  [t, yNew] = rk4_step(f1, t0, h, y);
end
