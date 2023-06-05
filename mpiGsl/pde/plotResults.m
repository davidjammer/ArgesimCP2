function plotResults(file)
% creates required plots

data = load(file);
t = data(2:end,1)';
y = data(2:end,2:end);
x = data(1,2:end);

par.L = x(end);
par.N = numel(x) - 1;
par.k = par.L / par.N ;
par.tEnd = t(end);
par.xP = [1/2, 3/4]*par.L;
par.tP = [0.5, 0.8]*par.tEnd;

subplot(2,1,1)
I = round(par.xP/par.k + 1);
plot(t, y(:,I))
grid("on")
title("excitation over time")
xlabel("t")
legend("x=0.25", "x=0.375", "Location", "best")

subplot(2,1,2)
u1 = y(t==par.tP(1),1:end);
u2 = y(t==par.tP(2),1:end);
plot(x, u1, x, u2)
grid("on")
title("excitation over space")
xlabel("x")
legend("t=5", "t=8", "Location", "best")

figure
m = size(y,1);
n = size(y,2);
[X,T] = meshgrid(x,t);
z = zeros(m,1);

surf(X, T, y, "EdgeColor", "none")
view(130,20)
colormap("jet")
hold("on")
contour3(X, T, y, "k-")
end
