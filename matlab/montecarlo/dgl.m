function [dydt] = dgl(t,y,k,d,m)

	dydt = zeros(2,1);
	dydt(1) = y(2);
	dydt(2) = -d / m * y(2) - k / m * y(1);
	
end

