%% Optimum design for -45dB return loss from pattern search optimization 
clear;
c0 = 3e8;
F = 5e9;
lamb = c0./F;
N = 5;
data = load('ga_V2L_ms3serv2_fl_V2.mat');
fmin = data.fmin2;

R = [fmin.r(1) fmin.r(1:N)];

% R = [0.02283 0.02283 0.06313 0.07772 0.09971 0.1467]; % From GA algo


% R = [0.02 0.02 0.02667 0.03333 0.04 0.04667 0.05333 0.06 0.06667 0.07333 0.09562 0.08667 0.09333 0.1009 0.1067 0.1156...
%     0.1491];
% N = length(R) - 1;
% Len = 0.3912;
% Len = fmin.r(N + 1);

Len = fmin.r(N+1);

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
title(['Optimum Horn, Aprture Radius = ', num2str(R(end)./lamb), ' \lambda ', 'Horn Lengt = ', num2str(Len./lamb), ' \lambda'], 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Horn Length cut', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Horn Radius cut', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% Iteration over the frequencies 

SP = [R(2:end) Len];

XMN_data = load('Xmn.mat');
str = XMN_data.Xmn;

fc_ = fc(R(2), 1, 1);

% <<<<<<< Updated upstream
% F = linspace(fc_(1)+fc_(1)./100, fc_(3), 20);
% =======
% F = linspace(fc_(1)+fc_(1)./100, 20e9, 50);
F = [4e9 6e9 8e9];
% >>>>>>> Stashed changes
% F = 5e9;
% objective = @(x) GSM_N_opt_allvar(SP(1:N), SP(N+1), x(1:length(F)), 0);
tic;
% parfor i = 1:length(F)

[Eth_, Eph_, Eco_, Exp_, CO_, XP_, E, th, ph, Max_Exp_diff] = MinXP_Goal_V2L_freq_fields(SP(1:N), SP(N+1), F, 20);



% save('Ga_SRR', 'SRR');

time_used = toc;
% end
%% plots

figure;
surf(F*1e-9, th(1, :)*180/pi, db(abs(squeeze(E(:, 1, :)).')), 'LineWidth', 2); shading flat;
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\theta (deg)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('E(\theta, 0) (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('E(\theta, 0) magnitude in dB', 'FontSize', 12, 'FontWeight', 'bold');
% RL1 = GSM_N_opt_allvar(SP(1:N), SP(N+1), 5e9, 20);


figure;
plot(F*1e-9, Max_Exp_diff, 'LineWidth', 2); shading flat;
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cross polar level (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Cross polar level from RL optimized antenna using GA', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% Directivity
zeta = 120*pi;
figure;
for i = 1:length(F)
    E_Ip = sqrt(abs(squeeze(Eth_(i, :, :))).^2 + abs(squeeze(Eph_(i, :, :))).^2);
    U_feedIp = abs(E_Ip).^2./(2 .* zeta);
    P_radIp_i = U_feedIp .* sin(th(1, :)) .*pi./180 .* pi./180;
    P_radIp = sum(sum(P_radIp_i));
    D_Ip(i, :, :) = 4 .* pi .* U_feedIp ./ P_radIp;
    hold on;
    plot(th(1, :)*180/pi, db(squeeze(D_Ip(i, 1, :)))./2, 'LineWidth', 2);
end

