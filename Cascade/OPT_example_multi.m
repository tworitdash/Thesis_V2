% clear;
% 
% % [STT, STR, SRT, SRR, N] = GSM_N(R, L, er, mur, F, k);
% c0 = 3e8;
% F = 5e9;
% 
% lamb = c0./F;
% 
% Len = 6 .* lamb;
% % 
% % R1 = 2e-2;
% % Rend = 2 .* lamb;
% 
% % funct = @(x) GSM_N_opt(R1, Rend, x(1), Len, F, 0);
% 
% 
% problem.objective =  @(x) GSM_N_opt(x(2), x(3), x(1), Len, F, 0);
% problem.x0 = [7 2e-2 4e-2];
% problem.lb = [7 2e-2 4e-2];
% problem.ub = [round(Len./(lamb./10)) 4e-2 8e-2];
% problem.solver = 'fmincon';
% 
% tic;
% 
% % problem.options = optimoptions(@ga,'PlotFcn', {'gaplotbestf', 'gaplotscores'},...
% %     'InitialPopulationMatrix', [round(Len./(lamb./40))], 'UseParallel', true);
% 
% problem.options = optimoptions(@fmincon,'PlotFcn', ...
%     {'optimplotx', 'optimplotfirstorderopt'},...
%      'UseParallel', true);
% 
% [x, fval, exf, ouput] = fmincon(problem);
% 
% time_opt = toc;
% 
% fmin.x = x;
% fmin.fval = fval;
% fmin.exf = exf;
% fmin.output = ouput;
% fmin.time_consumed = time_opt;
% 
% save('fmincon_3var_r1rendN_RL.mat', 'fmin');

%% Optimizing the radii of each section 


N = round(fmin.x(1));
R1 = fmin.x(2);
Rend = fmin.x(3);

R_test = linspace(R1, Rend, N);

problem2.objective = @(x) GSM_N_opt_radii([x(1:N)], Len, F, 0);
problem2.x0 = [linspace(R1, Rend, N)];

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
problem2.solver = 'fmincon';


tic;

% problem.options = optimoptions(@ga,'PlotFcn', {'gaplotbestf', 'gaplotscores'},...
%     'InitialPopulationMatrix', [round(Len./(lamb./40))], 'UseParallel', true);

problem2.options = optimoptions(@fmincon,'PlotFcn', ...
    {'optimplotx', 'optimplotfirstorderopt'},...
     'UseParallel', true);

[r, fval2, exf2, ouput2] = fmincon(problem2);

time_opt2 = toc;

fmin2.r = r;
fmin2.fval = fval2;
fmin2.exf = exf2;
fmin2.output = ouput2;
fmin2.time_consumed = time_opt2;

save('fmincon_Nvar_Rvec.mat', 'fmin2');

