% clear;
c0 = 3e8;
N = 5;
F = 3e9;

lamb = c0/F;

Len = lamb;

str = load('Xmn_azimuthal_inc_TE.mat');
str = str.xmn_TE;

xmn = str(1).xmn;

R1 = (xmn .* c0)./(2 .* pi .* F) + 1e-3;

% R1 = 2e-2;
Rend = 0.6579 .* lamb;
% Rend = 20e-2;

R_test = linspace(R1, Rend, N);


focal_length = 1; % Focal length of the reflector
theta_0 = 50 .* pi./180; % Subtended angle given from SKA
d = 4 .* focal_length .* tan(theta_0./2); % Reflectoor Diameter

% problem2.fitnessfcn = @(x) MinXP_Goal([x(1:N)], [l1 x(N+1:2*N-1)], F, 0, 10, 6);

problem2.fitnessfcn = @(x) Ga_opt_aper_eff([x(1:N)], x(N + 1), F, 5, focal_length, d);
problem2.nvars = N + 1;
% problem2.x0 = [linspace(R1, Rend, N)];

lb = zeros(N+1, 1);
ub = zeros(N+1, 1);

lb(1) = R1;
ub(1) = R1 + lamb/4;

for i =  2:N
    lb(i) = ub(i - 1) - lamb/2;
    ub(i) = R_test(i) + lamb/2;
end

lb(N + 1) = 0.5*lamb;
ub(N + 1) = 2*lamb;

problem2.lb = [lb.'];
problem2.ub = [ub.'];

IP = [R_test Len];

problem2.solver = 'ga';


tic;
% N2 = N+1:2*N-1;
% 
% for i = 1:length(N2)
%     L(i) = sum(r(N+1:N2(i)));
% end
% L = [l1 L];

% problem.options = optimoptions(@ga,'PlotFcn', {'gaplotbestf', 'gaplotscores'},...
%     'InitialPopulationMatrix', [round(Len./(lamb./40))], 'UseParallel', true);

% problem2.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
%     'InitialPopulationMatrix', [IP], 'UseParallel',...
%     true, 'FitnessLimit',  10^(-45/20));

problem2.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
    'InitialPopulationMatrix', [IP], 'UseParallel',...
    true, 'Fitnesslimit', -0.75);

[r, fval2, exf2, ouput2] = ga(problem2);

time_opt2 = toc;

fmin2.r = r;
fmin2.fval = fval2;
fmin2.exf = exf2;
fmin2.output = ouput2;
fmin2.time_consumed = time_opt2;

% save('ga_allvar_ms3serv2_wo_fl.mat', 'fmin2'); %wo_fl is for without fitness limit

save('ga_allvar_ms3serv2_maxeap.mat', 'fmin2'); %minxp is for minimum cross polarization
