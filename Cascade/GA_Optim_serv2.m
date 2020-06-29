clear;
c0 = 3e8;
F = 5e9;

lamb = c0./F;

Len = 6 .* lamb;

% R1 = 2e-2;
% Rend = 2 .* lamb;
fun = @(x) GSM_N_opt(x(2), x(3), x(1), Len, F, 0);

problem.fitnessfcn = fun;
problem.nvars = 3;
problem.lb = [7 2e-2 4e-2];
problem.ub = [round(Len./(lamb./10)) 4e-2 8e-2];
problem.IntCon = 1;
problem.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
    'InitialPopulationMatrix', [7 2e-2 5e-2], 'UseParallel',...
    true, 'FitnessLimit',  10^(-35/20));
tic;

[x, fval, exf, out] = ga(problem);

time_opt = toc;

fmin.x = x;
fmin.fval = fval;
fmin.exf = exf;
fmin.output = out;
fmin.time_consumed = time_opt;

save('ga_3var_r1rendN_RL_ms3serv2.mat', 'fmin');

%% 
N = round(fmin.x(1));
R1 = fmin.x(2);
Rend = fmin.x(3);

R_test = linspace(R1, Rend, N);
problem2.fitnessfcn = @(x) GSM_N_opt_radii([x(1:N)], Len, F, 0);
problem2.nvars = N;
% problem2.x0 = [linspace(R1, Rend, N)];

lb = zeros(N, 1);
ub = zeros(N, 1);

lb(1) = R1 - R1/5;
ub(1) = R1 + R1/5;

for i =  2:N
    lb(i) = ub(i - 1) - R_test(1)/5;
    ub(i) = R_test(i) + R_test(1)/5;
end

problem2.lb = [lb];
problem2.ub = [ub];


tic;

% problem.options = optimoptions(@ga,'PlotFcn', {'gaplotbestf', 'gaplotscores'},...
%     'InitialPopulationMatrix', [round(Len./(lamb./40))], 'UseParallel', true);

problem2.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
    'InitialPopulationMatrix', [[linspace(R1, Rend, N)]], 'UseParallel',...
    true);

[r, fval2, exf2, ouput2] = ga(problem2);

time_opt2 = toc;

fmin2.r = r;
fmin2.fval = fval2;
fmin2.exf = exf2;
fmin2.output = ouput2;
fmin2.time_consumed = time_opt2;

save('ga_Nvar_Rvec_ms3serv2.mat', 'fmin2');


