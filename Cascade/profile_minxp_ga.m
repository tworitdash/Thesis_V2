%% Optimum design for -45dB return loss from pattern search optimization 
clear;
c0 = 3e8;
F = 5e9;
lamb = c0./F;
N = 5;
data = load('ga_allvar_minxp_V2L.mat');
fmin = data.fmin2;

R = [fmin.r(1) fmin.r(1:N)];

% R = [0.02283 0.02283 0.06313 0.07772 0.09971 0.1467]; % From GA algo


% R = [0.02 0.02 0.02667 0.03333 0.04 0.04667 0.05333 0.06 0.06667 0.07333 0.09562 0.08667 0.09333 0.1009 0.1067 0.1156...
%     0.1491];
% N = length(R) - 1;
% Len = 0.3912;
Len = fmin.r(N + 1);

% Len = 0.3273;

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
F = linspace(fc_(1)+fc_(1)./100, fc_(3), 20);
% =======
% F = linspace(fc_(1)+fc_(1)./100, 20e9, 20);
% >>>>>>> Stashed changes
% F = 5e9;
% objective = @(x) GSM_N_opt_allvar(SP(1:N), SP(N+1), x(1:length(F)), 0);
tic;
% parfor i = 1:length(F)
% XP_level = MinXP_Goal_V2L_freq(SP(1:N), SP(N+1), F, 20);
F_ = [F(1) F(4) F(10) F(20)];
[Eth_, Eph_, Eco_, Exp_, CO_, XP_, E_, th, ph, Max_Exp_diff] = MinXP_Goal_V2L_freq_fields(SP(1:N), SP(N + 1), 2.6129e9, 10);
% figure;
% plot(F_*1e-9, Max_Exp_diff, 'LineWidth', 2);grid on;
% xlabel('Frequency (GHz)', 'FontSize', 16, 'FontWeight', 'bold');
% ylabel('Cross Polar level (dB)', 'FontSize', 16, 'FontWeight', 'bold');
% title('Cross polar level ga algo', 'FontSize', 16, 'FontWeight', 'bold');

% save('Max_XP', 'Max_Exp_diff');

%% Plots

i = 1;

Eco = squeeze(Eco_(i, :, :));
Exp = squeeze(Exp_(i, :, :));
E = squeeze(E_(i, :, :));

CO = squeeze(CO_(i, :, :));

XP = squeeze(XP_(i, :, :));

time_used = toc;
% end

% % RL1 = GSM_N_opt_allvar(SP(1:N), SP(N+1), 5e9, 20);
% end


figure(80);

hold on;
plot(th(1, :)*(180/pi), db((abs(E(1, :))/max(abs(E(1, :))))), '-.', 'LineWidth', 2);

hold on;
plot(th(90, :)*(180/pi), db((abs(E(90, :))/max(abs(E(90, :))))), '-.', 'LineWidth', 2);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs} (dB) normalized ', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\phi = 0', '\p
grid on;
ylim([-40 5]);

figure(82);
hold on;
plot(th(1, :)*(180/pi), db((abs(Eco(1, :))./max(max(abs(Eco))))), '-.', 'LineWidth', 2);
hold on;
plot(th(90, :)*(180/pi), db((abs(Exp(45, :))./max(max(abs(Eco))))), '-.', 'LineWidth', 2);

xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Far fields E_{co} and E_{xp} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Co and Cross polar fields', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
ylim([-40 5]);

% 
% figure(72);
% hold on;
% plot(th(1, :)*(180/pi), db((abs(CO(1, :))./max(max(abs(CO))))), '-.', 'LineWidth', 2);
% hold on;
% plot(th(90, :)*(180/pi), db((abs(XP(45, :))./max(max(abs(CO))))), '-.', 'LineWidth', 2);
% 
% xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Far fields E_{co} and E_{xp} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Co and Cross polar fields', 'FontSize', 12, 'FontWeight', 'bold');
% grid on;
% ylim([-40 5]);