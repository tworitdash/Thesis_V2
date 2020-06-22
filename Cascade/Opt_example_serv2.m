clear;

% [STT, STR, SRT, SRR, N] = GSM_N(R, L, er, mur, F, k);
c0 = 3e8;
F = 5e9;

lamb = c0./F;

Len = 6 .* lamb;

R1 = 2e-2;
Rend = 2 .* lamb;

fun = @(x) GSM_N_opt(R1, Rend, x(1), Len, F, 0);

% [x, fval] = ga(fun, 2, [], [], [], [], [2e-2 4e-2], [4e-2 8e-2]);

x0 = [2e-2, 4e-2, round(Len./(lamb./10))];

% a = fun(x0)

A = [];
b = [];
Aeq = [];
beq = [];
% lb = [2e-2, 4e-2, 7];
% ub = [4e-2, 8e-2, round(Len./(lamb./40))];

lb = [7];
ub = [round(Len./(lamb./40))];

tic;

%[x, fval, ef, output, lambda] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub);

%options = optimoptions('simulannealbnd','Display', 'iter');

options = optimoptions(@ga,'PlotFcn', {'gaplotbestf', 'gaplotscores'},...
    'InitialPopulationMatrix', [round(Len./(lamb./40))], 'UseParallel', true, 'FitnessLimit',  10^(-15/20));

[x, fval, exf, ouput, population] = ga(fun, 1, A, b, Aeq, beq, lb, ub, [], options);

time_opt = toc;


