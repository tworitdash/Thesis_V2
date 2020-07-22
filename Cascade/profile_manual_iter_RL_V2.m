%% Optimum design for -45dB return loss from pattern search optimization 
clear;
c0 = 3e8;
F = 5e9;
lamb = c0./F;
N = 16;
data = load('fmincon_V2L_ms3serv2_fl.mat');
fmin = data.fmin2;

R = [fmin.r(1) fmin.r(1:N).'];

% R = [0.02 0.02 0.02667 0.03333 0.04 0.04667 0.05333 0.06 0.06667 0.07333 0.09562 0.08667 0.09333 0.1009 0.1067 0.1156...
%     0.1491];
% N = length(R) - 1;
% Len = 0.3912;
Len = fmin.r(N + 1);

l1 = lamb./4;
L = [l1 ones(1, N - 1)./N .* Len];

for i = 1:N
    L_axis(i) = sum(L(1:i));
end
L_axis = [0 L_axis];
figure;
plot(L_axis, R, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] );
hold on;
plot(L_axis, -R, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] );

grid on;

%% Iteration over the frequencies 

SP = [R(2:end) Len];

XMN_data = load('Xmn.mat');
str = XMN_data.Xmn;

fc_ = fc(R(2), 1, 1);

F = linspace(fc_(1)+fc_(1)./100, fc_(3), 1000);
% F = 5e9;
% objective = @(x) GSM_N_opt_allvar(SP(1:N), SP(N+1), x(1:length(F)), 0);

% parfor i = 1:length(F)
 RL = GSM_N_opt_allvar(SP(1:N), SP(N+1), F, 20);
% end

figure;
plot(F*1e-9, RL, 'LineWidth', 2);grid on;
xlabel('Frequency (GHz)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Return loss (dB)', 'FontSize', 16, 'FontWeight', 'bold');
title('RL of patternsearch algo', 'FontSize', 16, 'FontWeight', 'bold');
% RL1 = GSM_N_opt_allvar(SP(1:N), SP(N+1), 5e9, 20);
% end
