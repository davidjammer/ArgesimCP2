function dydt = ode_matrix(t,y,par, A)
% DGL of the discretized PDE, using matrix-vector style
% y is Composite, column vector of size 2*lsize, [u; ut]
%      contains left/right boundary value on first/last task
  RIGHT = 1;
  LEFT = 2;
  n = numel(y)/2;             % == lsize

  if numlabs == 1
    u = y(1:n);
    ut = y(n+1:end);
    dx = ut;
    dv = A*u;
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
    dv = [0; A(2:n,:)*u];
  elseif labindex == numlabs
    dx = [ut(1:end-1); 0];
    dv = [A(2:n,:)*u; 0];
  else
    dx = ut;    
    dv = A(2:n+1,:)*u;
  end
  
  dydt = [dx; dv];
end
