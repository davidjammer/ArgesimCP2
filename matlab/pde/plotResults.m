function plotResults(y, par)
% creates required plots
t = 0:par.dtOut:par.tEnd;
n = par.N - 1;

subplot(2,1,1)
I = round(par.xP/par.k + 1);
plot(t, y(:,I))
grid("on")
title("excitation over time")
xlabel("t")
legend("x=0.25", "x=0.375", "Location", "best")

subplot(2,1,2)
x = 0:par.k:par.L;
u1 = [0, y(t==par.tP(1), 1:n), 0];
u2 = [0, y(t==par.tP(2), 1:n), 0];
plot(x, u1, x, u2)
grid("on")
title("excitation over space")
xlabel("x")
legend("t=5", "t=8", "Location", "best")

figure
m = size(y,1);
[X,T] = meshgrid(x,t');
z = zeros(m,1);
Y = [z, y(:,1:n), z];
surf(X, T, Y, "EdgeColor", "none")
view(130,20)
colormap("jet")
hold("on")
contour3(X, T, Y, "k-")
end
