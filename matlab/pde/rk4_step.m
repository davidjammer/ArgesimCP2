function [t, yNew] = rk4_step(f, t0, h, y)
	t = t0;
	k1 = h*f(t, y);
	k2 = h*f(t + h/2, y + k1/2);
	k3 = h*f(t + h/2, y + k2/2);
	k4 = h*f(t + h, y + k3);

	yNew = y + (k1 + 2*k2 + 2*k3 + k4)/6;
	t = t + h;
end
