function tE = run(location, type, n, dt, N)
% run pde function on parallel cluster
% parameters:
%   location    "local" | "cluster"
%   type        "loop" | "matrix" | "dist" | "dist1" | "distSlow" 
%   n           number of tasks
%   dt          time step size of output values
%   N           number of intervals for space discretisation
  if nargin < 4
    dt = 0.01;
    N = 500;
  end
 
  if location == "local"
    fprintf("run local with %d Worker(s)\n", n);
    pool = parpool(n);
    switch type
      case "loop"     % explicit loops
        [u, par, tE] = pde_loop(dt, N);
      case "matrix"   % matrix-vector style
        [u, par, tE] = pde_matrix(dt, N);
      case "dist"     % distributed matrices
        [u, par, tE] = pde_dist(dt, N);
      case "dist1"     % distributed matrices, tests
        [u, par, tE] = pde_dist1(dt, N);
      otherwise       % distributed matrices, slow version
        [u, par, tE] = pde_distslow(dt, N);
    end
    delete(pool);
  elseif location == "cluster"
    fprintf("run on cluster with %d Worker(s)\n", n);
    cluster = parcluster("Seneca");
    switch type
      case "loop"     % explicit loops
        job = batch(cluster, @pde_loop, 3, {dt, N}, "Pool", n, ...
          "AttachedFiles", ...
          ["rk4_step", "ode_loop", "setParams", "setInitialCondition"], ...
          "AutoAddClientPath", false);
      case "matrix"   % matrix-vector style
        job = batch(cluster, @pde_matrix, 3, {dt, N}, "Pool", n, ...
          "AttachedFiles", ...
          ["rk4_step", "ode_matrix", "setParams", "setInitialCondition"], ...
          "AutoAddClientPath", false);
      case "dist"     % distributed matrices
        job = batch(cluster, @pde_dist, 3, {dt, N}, "Pool", n, ...
          "AttachedFiles", ["setParams", "setInitialCondition"], ...
          "AutoAddClientPath", false);
      case "dist1"     % distributed matrices, tests
        job = batch(cluster, @pde_dist1, 3, {dt, N}, "Pool", n, ...
          "AttachedFiles", ["rk4_step", "setParams", "setInitialCondition"], ...
          "AutoAddClientPath", false);
      otherwise       % distributed matrices, slow version
        job = batch(cluster, @pde_distslow, 3, {dt, N}, "Pool", n, ...
          "AttachedFiles", ["rk4_step", "setParams", "setInitialCondition"], ...
          "AutoAddClientPath", false);
    end

    wait(job);
    results = fetchOutputs(job);
    u = results{1};
    par = results{2};
    tE = results{3};

    delete(job);
  end
  
  fprintf("Elapsed time is %f seconds\n", tE);
  plotResults(u, par);
end
