c0 = 3e8;
N = 5;
F = 5e9;

lamb = c0/F;

Len = 6 .* lamb;

R1 = 2e-2;
Rend = 2 .* lamb;

R_test = linspace(R1, Rend, N);

% l1 = lamb/4;
% l = ones(N-1, 1) * Len/N;

problem2.fitnessfcn = @(x) GSM_N_opt_allvar_V2([x(1:N)], x(N+1), F, 0);
problem2.nvars = N + 1;
% problem2.x0 = [linspace(R1, Rend, N)];

lb = zeros(N+1, 1);
ub = zeros(N+1, 1);

lb(1) = R1;
ub(1) = R1 + lamb./2;

for i =  2:N
    lb(i) = ub(i - 1) - lamb./2;
    ub(i) = R_test(i) + lamb./2;
end

lb(N + 1) = 5 .* lamb;
ub(N + 1) = 8  .* lamb;

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
    true, 'FitnessLimit', -25);

[r, fval2, exf2, ouput2] = ga(problem2);

time_opt2 = toc;

fmin2.r = r;
fmin2.fval = fval2;
fmin2.exf = exf2;
fmin2.output = ouput2;
fmin2.time_consumed = time_opt2;

save('ga_V2L_ms3serv2_fl_V2.mat', 'fmin2'); %wo_fl is for without fitness limit
