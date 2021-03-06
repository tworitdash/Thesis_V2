c0 = 3e8;
N = 16;
F = 5e9;

lamb = c0/F;

Len = 6 .* lamb;

R1 = 2e-2;
Rend = 2 .* lamb;

R_test = linspace(R1, Rend, N);

l1 = lamb/4;
l = ones(N-1, 1) * Len/N;

problem2.objective = @(x) GSM_N_opt_allvar([x(1:N)], [l1 x(N+1:2*N-1)], F, 0);
% problem2.nvars = 2*N - 1;
% problem2.x0 = [linspace(R1, Rend, N)];

lb = zeros(2.*N-1, 1);
ub = zeros(2.*N-1, 1);

lb(1) = R1;
ub(1) = R1 + R1/5;

for i =  2:N
    lb(i) = ub(i - 1) - R_test(1)/5;
    ub(i) = R_test(i) + R_test(1)/5;
end

for m = N+1:2*N-1
    lb(m) = l(m - N) - lamb/5;
    ub(m) = l(m - N) + lamb/5;
end

problem2.lb = [lb.'];
problem2.ub = [ub.'];

problem2.x0 = [R_test l.'];
problem2.solver = 'fmincon';


tic;
N2 = N+1:2*N-1;

for i = 1:length(N2)
    L(i) = sum(r(N+1:N2(i)));
end
L = [l1 L];

% problem.options = optimoptions(@ga,'PlotFcn', {'gaplotbestf', 'gaplotscores'},...
%     'InitialPopulationMatrix', [round(Len./(lamb./40))], 'UseParallel', true);

% problem2.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
%     'InitialPopulationMatrix', [IP], 'UseParallel',...
%     true, 'FitnessLimit',  10^(-45/20));

problem2.options = optimoptions(@fmincon,'PlotFcn', ...
    {'optimplotx', 'optimplotfirstorderopt'},...
     'UseParallel', true, 'MaxFunEval', Inf, 'MaxIter', Inf);

[r, fval2, exf2, ouput2] = fmincon(problem2);

time_opt2 = toc;

fmin2.r = r;
fmin2.fval = fval2;
fmin2.exf = exf2;
fmin2.output = ouput2;
fmin2.time_consumed = time_opt2;

save('fmincon_allvar_ms3serv2_V2.mat', 'fmin2');