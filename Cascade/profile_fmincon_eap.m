%% Optimum design for -45dB return loss from pattern search optimization 
clear;
c0 = 3e8;
F = 3e9;
lamb = c0./F;
N = 5;
data = load('fmincon_allvar_ms3serv2_maxeap.mat');
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
title(['Optimum Horn, Aprture Radius = ', num2str(R(end)./lamb), ' \lambda ', 'Horn Length = ', num2str(Len./lamb), ' \lambda'], 'FontSize', 12, 'FontWeight', 'bold');
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
F = linspace(fc_(1)+fc_(1)./100, fc_(11), 50);
% >>>>>>> Stashed changes
% F = 5e9;
% objective = @(x) GSM_N_opt_allvar(SP(1:N), SP(N+1), x(1:length(F)), 0);
tic;
% parfor i = 1:length(F)
[Aper_n] = Ga_opt_aper_eff(SP(1:N), SP(N+1), F, 5);


% save('Ga_SRR', 'SRR');
% save('Ga_RL', 'RL');

time_used = toc;
% end

figure;
plot(F*1e-9, -Aper_n*100, 'LineWidth', 2);grid on;
xlabel('Frequency (GHz)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Aperture efficiency (dB)', 'FontSize', 16, 'FontWeight', 'bold');
title('Aperture efficiency of ga algo', 'FontSize', 16, 'FontWeight', 'bold');