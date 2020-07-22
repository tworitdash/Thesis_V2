%% Optimum design for -45dB return loss from pattern search optimization 
clear;
c0 = 3e8;
F = 5e9;
lamb = c0./F;

N = 16;
data = load('ga_V2L_ms3serv2_fl.mat');
r = data.fmin2.r;
R = [r(1) r(1:N)];
% R = [0.02253 0.02253 0.03353 0.03917 0.05562 0.05946 0.05901 0.07181 0.08365 0.1009 0.08275 0.08619 0.1107 0.1097 0.1106 0.1297...
%     0.1182];
% N = length(R) - 1;
Len = r(N + 1);

l1 = lamb./4;
L = [l1 ones(1, N - 1)./N .* Len];

for i = 1:N
    L_axis(i) = sum(L(1:i));
end
L_axis = [0 L_axis];
figure(1);
hold on;
plot(L_axis, R, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] );
hold on;
plot(L_axis, -R, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] );

grid on;

%% Iteration over the frequencies 

SP = [R(2:end) Len];

XMN_data = load('Xmn.mat');
str = XMN_data.Xmn;

fc_ = fc(R(2), 1, 1);

F = linspace(fc_(1)+fc_(1)./100, fc_(2), 1000);
% F = 5e9;
% objective = @(x) GSM_N_opt_allvar(SP(1:N), SP(N+1), x(1:length(F)), 0);

% for i = 1:length(F)
RL_ga = GSM_N_opt_allvar(SP(1:N), SP(N+1), F, 20);
% end

figure(2);
hold on;
plot(F*1e-9, RL_ga, 'LineWidth', 2); grid on;
xlabel('Frequency (GHz)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Return loss (dB)', 'FontSize', 16, 'FontWeight', 'bold');
title('RL of Ga algo', 'FontSize', 16, 'FontWeight', 'bold');


