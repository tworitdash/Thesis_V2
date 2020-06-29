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

save('ga_3var_r1rendN_RL.mat', 'fmin');