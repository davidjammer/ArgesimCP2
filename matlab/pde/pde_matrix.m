function [uAll, p, tE] = pde_matrix(dt, N)
  % partial differential equation example
  % mimics MPI version, but uses matrix-vector operations
  %   dt          time step size of output values
  %   N           number of intervals for space discretisation
  if nargin < 2
    dt = 0.01;
    N = 500;
  end
  
  tic;
  p = setParams(dt, N);
  y0 = setInitialCondition(p);

  % matrix for ode
  fA = @(N) (p.nu*p.nu)/(p.k*p.k) * ...
      (diag(ones(N,1), -1) - 2*diag(ones(N+1, 1)) + diag(ones(N,1), 1));
  
  spmd
    codist = codistributor1d(2, codistributor1d.unsetPartition, [1,N+1]);
    [imin, imax] = globalIndices(codist, 2);
    lsize = imax - imin + 1;

    if numlabs == 1
      A = fA(lsize - 1);
    elseif labindex == 1 || labindex == numlabs
      A = fA(lsize);
    else
      A = fA(lsize + 1);      
    end
    
    yLoc = [y0(imin:imax); y0(imin+(N+1):imax+(N+1))];
    t = 0;
    uLoc = zeros(p.nSteps + 1, lsize);   % collects u values
    uLoc(1,:) = yLoc(1:lsize)';
    for step=1:p.nSteps
      labBarrier;
      [t, yLoc] = rk4_step1(@ode_matrix, t, dt, yLoc, p, A);
      uLoc(step+1,:) = yLoc(1:lsize)';
    end
  end

  uAll = zeros(p.nSteps + 1, p.N + 1);
  for I=1:length(uLoc)
    uAll(:,imin{I}:imax{I}) = uLoc{I};
  end

  tE = toc;
end

function [t, yNew] = rk4_step1(f, t0, h, y, p, A)
  f1 = @(t,y) f(t, y, p, A);
  [t, yNew] = rk4_step(f1, t0, h, y);
end
