function t = monte_seq(nReps, nSteps)
% Monte-Carlo example, scalar
% parameters:
%    nReps     number of repetitions with varying d
%    nSteps    number of time steps

if nargin < 2
  nReps = 1000;
  nSteps = 200;
end

tic;
% given parameters
d0 = 800;
d1 = 1200;
x0 = 0;
v0 = 0.1;
tEnd = 2;

% other parameters
dtOut = tEnd/nSteps;
tP = (0:dtOut:tEnd)';

% compute all d's
d = (d1 - d0)*rand(nReps, 1) + d0;

% compute all results and the mean
xM = zeros(length(tP), 1);
for I=1:nReps
  f = @(t,y) dampOde(t,y,d(I));
  [tP,y] = ode45(f, tP, [x0, v0]);
  xM = xM + y(:,1);
end
xM = xM/nReps;
t = toc;

% plot mean result
plot(tP, xM)
grid("on")

%-------------------------------------------------------------
function dy = dampOde(t, y, d)
% DGL of the damped oscillator with damping constant d
m = 450;
k = 9000;

dy1 = y(2);
dy2 = -(k/m)*y(1) - (d/m)*y(2);

dy = [dy1; dy2];
