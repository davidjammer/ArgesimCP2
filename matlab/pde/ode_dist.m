function dydt = ode_para3(t,y,par)
% DGL of the discretized PDE, using distributed matrices
% y is Composite, column vector of size 2*lsize, [u; ut]
%      contains left/right boundary value on first/last task
  n = numel(y)/2;             % == lsize
  nuk = (par.nu*par.nu)/(par.k*par.k);
  fA = @(N) diag(ones(N,1), -1) - 2*diag(ones(N+1, 1)) + diag(ones(N,1), 1);

  if numlabs == 1
    A = fA(n-1);
    u = y(1:n);
    ut = y(n+1:end);
    dx = ut;
    dv = nuk*A*u;
    dydt = [dx; dv];
    return;
  end

  % create local variables including ghost points (only for u)
  if labindex == 1
    u = [y(1:n);0];
    ut = y(n+1:end);
  elseif labindex == numlabs
    u = [0;y(1:n)];
    ut = y(n+1:end);
  else
    u = [0;y(1:n);0];
    ut = y(n+1:end);
  end

  labBarrier;
  RIGHT = 1;
  LEFT = 2;
  if labindex == 1
    labSend(u(n), labindex+1, RIGHT);
    u(n+1) = labReceive(labindex+1, LEFT);
  elseif labindex == numlabs
    labSend(u(2), labindex-1, LEFT);
    u(1) = labReceive(labindex-1, RIGHT);
  else
    labSend(u(n+1), labindex+1, RIGHT);
    labSend(u(2), labindex-1, LEFT);
    u(1) = labReceive(labindex-1, RIGHT);
    u(n+2) = labReceive(labindex+1, LEFT);
  end

  if labindex == 1
    dx = [0; ut(2:end)];
    A = fA(n);
    dv = [0; nuk*A(2:n,:)*u];
  elseif labindex == numlabs
    dx = [ut(1:end-1); 0];
    A = fA(n);
    dv = [nuk*A(2:n,:)*u; 0];
  else
    dx = ut;    
    A = fA(n+1);
    dv = nuk*A(2:n+1,:)*u;
  end
  
  dydt = [dx; dv];
end
