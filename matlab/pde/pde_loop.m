function [uAll, p, tE] = pde_loop(dt, N)
  % partial differential equation example
  % mimics MPI version
  %   dt          time step size of output values
  %   N           number of intervals for space discretisation
  if nargin < 2
    dt = 0.01;
    N = 500;
  end
  
  tic;
  p = setParams(dt, N);
  y0 = setInitialCondition(p);
  fOde = @(t,y) ode_loop(t,y,p);

  spmd
    codist = codistributor1d(2, codistributor1d.unsetPartition, [1,N+1]);
    [imin, imax] = globalIndices(codist, 2);
    lsize = imax - imin + 1;
        
    yLoc = [y0(imin:imax); y0(imin+(N+1):imax+(N+1))];
    t = 0;
    uLoc = zeros(p.nSteps + 1, lsize);   % collects u values
    uLoc(1,:) = yLoc(1:lsize)';
    for step=1:p.nSteps
      labBarrier;
      [t, yLoc] = rk4_step(fOde, t, dt, yLoc);
      uLoc(step+1,:) = yLoc(1:lsize)';
    end
  end

  uAll = zeros(p.nSteps + 1, p.N + 1);
  for I=1:length(uLoc)
    uAll(:,imin{I}:imax{I}) = uLoc{I};
  end

  tE = toc;
end
