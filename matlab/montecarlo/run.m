function t = run(location, type, n, nReps, nSteps)
% run montecarlo function on parallel cluster
% parameters:
%   location    "local" | "cluster"
%   type        "parfor" | "spmd" | "parsim"
%   n           number of tasks
%   nReps       number of repetitions (d values)
%   nSteps      number of time steps per computation

  if nargin < 5
    nReps = 1000;
    nSteps = 200;
  end

  if location == "local"
    fprintf("run local with %d Worker(s)\n", n);
    pool = parpool(n);
    if type == "parsim"
      [tvec, ymean, t] = monte_parsim(nReps, nSteps);
    elseif type == "parfor"
      [tvec, ymean, t] = monte_parfor(nReps, nSteps);
    else
      [tvec, ymean, t] = monte_spmd(nReps, nSteps);
    end
    delete(pool);
  elseif location == "cluster"
    fprintf("run on cluster with %d Worker(s)\n", n);
    cluster = parcluster("Seneca");
    if type == "parsim"
      job = batch(cluster, @monte_parsim, 3, {nReps, nSteps}, ...
        "Pool", n, "AttachedFiles", "ddm", "AutoAddClientPath", false);
    elseif type == "parfor"
      job = batch(cluster, @monte_parfor, 3, {nReps, nSteps}, ...
        "Pool", n, "AttachedFiles", "dgl", "AutoAddClientPath", false);
    else
      job = batch(cluster, @monte_spmd, 3, {nReps, nSteps}, ...
        "Pool", n, "AttachedFiles", "dgl", "AutoAddClientPath", false);
    end

    wait(job);
    results = fetchOutputs(job);
    tvec = results{1};
    ymean = results{2};
    t = results{3};

    delete(job);
  end
  
  fprintf("Elapsed time is %f seconds\n", t);
  plot(tvec, ymean, "r-", [0, tvec(end)], [0, 0], "k-")
end
