clear;

% [STT, STR, SRT, SRR, N] = GSM_N(R, L, er, mur, F, k);
c0 = 3e8;
F = 5e9;

lamb = c0./F;

Len = 6 .* lamb;
% 
% R1 = 2e-2;
% Rend = 2 .* lamb;

% funct = @(x) GSM_N_opt(R1, Rend, x(1), Len, F, 0);


problem.objective =  @(x) GSM_N_opt(x(2), x(3), x(1), Len, F, 0);
problem.x0 = [round(Len./(lamb./40)) 2e-2 4e-2];
problem.lb = [7 2e-2 4e-2];
problem.ub = [round(Len./(lamb./40)) 4e-2 8e-2];
problem.solver = 'fmincon';

tic;

% problem.options = optimoptions(@ga,'PlotFcn', {'gaplotbestf', 'gaplotscores'},...
%     'InitialPopulationMatrix', [round(Len./(lamb./40))], 'UseParallel', true);

problem.options = optimoptions(@fmincon,'PlotFcn', ...
    {'optimplotx', 'optimplotfirstorderopt'},...
     'UseParallel', true);

[x, fval, exf, ouput] = fmincon(problem);

time_opt = toc;


