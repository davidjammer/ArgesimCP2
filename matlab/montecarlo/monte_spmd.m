function [tvec, ymean, t] = monte_spmd(nReps, nSteps)
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

  spmd
    l_rep = fix(nReps / numlabs);
    if mod(nReps,numlabs) > 0 && (numlabs-mod(nReps,numlabs) < labindex)
      l_rep = l_rep+1;
    end
    
    D = rand(1,l_rep) * (D_MAX - D_MIN) + D_MIN;
    
    y = zeros(numel(tvec),1);
    for i=1:l_rep
      [~, y_temp] = dglcall(tvec, y0,K,D(i),M);
      y = y + y_temp(:,1);
    end
    ysum = gop(@plus,y,1);
  end

  ymean = ysum{1}./nReps;
  t=toc;
end

function [t,y] = dglcall(tvec,y0,K,D,M)
  [t, y] = ode45(@(t,y) dgl(t,y,K,D,M), tvec, y0);
end