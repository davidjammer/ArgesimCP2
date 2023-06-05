function y0 = setInitialCondition(p)
  % returns inital condition of state vector + edge points
  n = p.N - 1;         % number of unknowns
  y0 = zeros(2*n + 4, 1);  % u and u_t, incl. edge points

  iM = (n+1)/2;        % index of middle point
  I1 = (2:iM+1)'; 
  y0(I1) = 2*(p.h/p.L)*p.k*I1;
  I2 = ((iM+2):n+1)'; 
  y0(I2) = 2*p.h*(1 - (p.k/p.L)*I2);
end
