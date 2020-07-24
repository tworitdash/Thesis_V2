c0 = 3e8;

data = load('ga_WKB_minxp.mat', 'fmin2'); 

F = 5e9;
lamb = c0./F;

data = data.fmin2;

Rbase = data.r(3);
Len = data.r(2);
Slope = data.r(1);

R_axis = [Rbase Rbase Rbase + Slope * Len];
L_axis = [0 lamb/4 lamb/4 + Len];

figure;
plot(L_axis, R_axis, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] ); hold on; grid on; plot(L_axis, -R_axis, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840]);
title(['Optimum Horn, Aprture Radius = ', num2str(R_axis(end)./lamb), ' \lambda ', 'Horn Length = ', num2str(Len./lamb), ' \lambda'], 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Horn Length cut', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Horn Radius cut', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% 
SP = [R_axis(2:end) Len];

XMN_data = load('Xmn.mat');
str = XMN_data.Xmn;

fc_ = fc(R_axis(2), 1, 1);

Freq = linspace(fc_(1)+fc_(1)./100, fc_(3), 1000);
% F = 5e9;
% objective = @(x) GSM_N_opt_allvar(SP(1:N), SP(N+1), x(1:length(F)), 0);

parfor i = 1:length(F)
 Max_Exp_diff(i) = MinXP_Goal_V2L(SP(1:2), Len, 5e9, Freq(i), Freq(end));
end

figure;
plot(F*1e-9, Max_Exp_diff, 'LineWidth', 2);grid on;
xlabel('Frequency (GHz)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('RL in (dB)', 'FontSize', 16, 'FontWeight', 'bold');
title('WKB + Patternsearch algo for MinXP', 'FontSize', 16, 'FontWeight', 'bold');