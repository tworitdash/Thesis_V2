c0 = 3e8;
N = 40;
F = 5e9;

lamb = c0/F;

Len = 6 .* lamb;

R1 = 2e-2;
Rend = 2 .* lamb;

R_test = linspace(R1, Rend, N);


problem2.objective = @(x) MinXP_Goal_V2L([x(1:N)], x(N+1), F, 0);
% problem2.nvars = N + 1;
% problem2.x0 = [linspace(R1, Rend, N)];

lb = zeros(N+1, 1);
ub = zeros(N+1, 1);

lb(1) = R1;
ub(1) = R1 + R1/5;

for i =  2:N
    lb(i) = ub(i - 1) - lamb/2;
    ub(i) = R_test(i) + lamb/2;
end
ub(N + 1) = 10 .* lamb;
lb(N + 1) = 5 .* lamb; 

% for m = N+1:2*N-1
%     lb(m) = l(m - N) - lamb/5;
%     ub(m) = l(m - N) + lamb/5;
% end

problem2.lb = [lb.'];
problem2.ub = [ub.'];

IP = [R_test Len];
problem2.x0 = IP;

problem2.solver = 'fmincon';


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

% problem2.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
%     'InitialPopulationMatrix', [IP], 'UseParallel',...
%     true, 'Fitnesslimit', -20);

problem2.options = optimoptions(@fmincon,'PlotFcn', ...
    {'optimplotx', 'optimplotfirstorderopt'},...
     'UseParallel', true, 'MaxIterations', 200);

[r, fval2, exf2, ouput2] = fmincon(problem2);

time_opt2 = toc;

fmin2.r = r;
fmin2.fval = fval2;
fmin2.exf = exf2;
fmin2.output = ouput2;
fmin2.time_consumed = time_opt2;

% save('ga_allvar_ms3serv2_wo_fl.mat', 'fmin2'); %wo_fl is for without fitness limit

save('fmincon_allvar_ms3serv2_minxp_V2L_serv2.mat', 'fmin2'); %minxp is for minimum cross polarization
