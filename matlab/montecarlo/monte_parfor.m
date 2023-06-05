function [tvec, ymean, t] = monte_parfor(nReps, nSteps)
  tic;
  TSTART = 0;
  TEND = 2;

  K = 9000;
  M = 450;
  D_MIN = 800;
  D_MAX = 1200;

  h = TEND/nSteps;
  tvec = (TSTART:h:TEND)';
  y0 = [0, 0.1];

  D = rand(1,nReps) * (D_MAX - D_MIN) + D_MIN;
  ysum = zeros(size(tvec));
  parfor i=1:numel(D)
    [~, yout] = ode45(@(t,y) dgl(t,y,K,D(i),M), tvec, y0);
    ysum = ysum + yout(:,1);
  end

  ymean = ysum./nReps;
  t=toc;
