clear;
c0 = 3e8;
F = 5e9;

lamb = c0./F;

Len = 6 .* lamb;

R1 = 2e-2;
Rend = 2 .* lamb;

fun = @(x) GSM_N_opt(R1, Rend, x(1), Len, F, 0);

problem.fitnessfcn = fun;
problem.nvars = 1;
problem.lb = [7];
problem.ub = [round(Len./(lamb./40))];
problem.IntCon = 1;
problem.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
    'InitialPopulationMatrix', [round(Len./(lamb./40))], 'UseParallel',...
    true, 'FitnessLimit',  10^(-15/20));
tic;

[x, fval, exf, out] = ga(problem);

time_opt = toc;


