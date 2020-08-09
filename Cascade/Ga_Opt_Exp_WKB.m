clear;
c0 = 3e8;

F = 5e9;
lamb = c0/F;

rr = 2e-2;

er = 1; mur = 1;

problem.fitnessfcn = @(x)  WKB_Exp(x(1), x(2), x(3), F, er, mur);
problem.solver = 'ga';
problem.nvars = 3;

lb(3) = rr;
ub(3) = rr + lamb;

lb(2) = 4 .* lamb;
ub(2) = 10 .* lamb;

lb(1) = (lamb - lb(3))./lb(2);
ub(1) = (4 .* lamb - ub(3))./ub(2);

problem.lb = lb;
problem.ub = ub;

IP = [(2 .* lamb - lb(3))./(6 .* lamb) 6.*lamb rr];


problem.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
    'InitialPopulationMatrix', IP, 'UseParallel',...
    true, 'FitnessLimit', -25);
tic;

[r, fval2, exf2, ouput2] = ga(problem);

time_opt2 = toc;

fmin2.r = r;
fmin2.fval = fval2;
fmin2.exf = exf2;
fmin2.output = ouput2;
fmin2.time_consumed = time_opt2;

save('ga_WKB_minxp.mat', 'fmin2'); 