function par = setParams(dt, N)
  % returns struct of parameters for pde task
  if nargin < 2
    dt = 0.01;
    N = 500;
  end
  
  par.nu = 0.06;            % 0.6 in [1] is wrong
  par.L = 0.5;
  par.h = 0.05;
  par.N = N;

  par.tEnd = 10;
  par.dtOut = dt;
  par.xP = [1/2, 3/4]*par.L;
  par.tP = [0.5, 0.8]*par.tEnd;

  par.k = par.L/par.N;
  par.nSteps = fix(par.tEnd/par.dtOut);
end
