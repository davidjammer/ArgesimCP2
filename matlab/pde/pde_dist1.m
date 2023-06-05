function [uAll, p, tE] = pde_dist1(dt, N)
  % partial differential equation example
  % very slow, distribution is wrong! Looks like Matlab bug!!
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
    
    % distribution of u, ut
    % example values for:  run("local", "dist", 2, 0.1, 20)
    % V1: distributes [19+0], [0+19]
    %u = codistributed(u0, codistributor1d(1));
    %ut = codistributed(ut0, codistributor1d(1));
    % V2: distributes [19+0], [0+19]
    codimA = getCodistributor(A);
    codimU = codistributor1d(1, codimA.Partition);
    u = codistributed(u0, codimU);
    ut = codistributed(ut0, codimU);
    % V3: distributes [19+0], [0+19] !!
    %codimU = codistributor1d(1, [10,9]);
    %u = codistributed(u0, codimU);
    %ut = codistributed(ut0, codimU);
    
    uAll = zeros(p.nSteps + 1, N+1, "codistributed");  % collects u values
    uAll(1,2:N) = u0(1:N-1,:)';
    
    for step=1:p.nSteps
      [t, u, ut] = rk4_step1(@pdeOde, t, dt, u, ut, A);
      uAll(step+1,2:N) = u(1:N-1,:)';
    end
  end
  uAll = gather(uAll);
  tE = toc;

  spmd, cu = getCodistributor(u); cut = getCodistributor(ut); end;
  cu = cu{1}; cut = cut{1};
  fprintf("\nPartitioning of u: ");
  fprintf(" %d ", cu.Partition);
  fprintf("\nPartitioning of ut: ");
  fprintf(" %d ", cut.Partition);
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

function [t, u, ut] = rk4_step1(f, t0, h, u, ut, A) 
  f1 = @(t,y) f(t, y, A);
  [t, yNew] = rk4_step(f1, t0, h, [u; ut]);
  
  n = numel(u);
  u = yNew(1:n);
  ut = yNew(n+1:end);
end
