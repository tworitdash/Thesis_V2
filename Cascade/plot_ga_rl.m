MM = load('Ga_SRR.mat');
SRR = MM.SRR;

MM_RL = load('Ga_RL.mat');
RL = MM_RL.RL;

% F_feko = linspace(4e9, 8e9, 36);
F_feko = linspace(4e9, 8e9, 36);
% F_feko = F;

data = read(rfdata.data,'../../Thesis/3wg/gsm/RL_Ga.s2p');
s_params = extract(data,'S_PARAMETERS');

RL_feko = db(abs(squeeze(s_params(1, 1, :)) + squeeze(s_params(1, 2, :))).^2)./2;

figure;
hold on;

plot(F * 1e-9, RL, 'LineWidth', 2); hold on;
plot(F_feko * 1e-9, RL_feko, 'LineWidth', 2); hold on;
grid on;
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Return loss (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('RL of fmincon algo', 'FontSize', 12, 'FontWeight', 'bold');
legend({'MM', 'FEKO'}, 'FontSize', 12, 'FontWeight', 'bold');

figure;

plot(F * 1e-9, db(abs(squeeze(SRR(:, 9, 1)))), 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, db(abs(squeeze(s_params(2, 1, :)))), 'LineWidth', 2); grid on;
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S11 of TE_{11} with TE_{12} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('S11 in dB', 'FontSize', 12, 'FontWeight', 'bold');
legend({'MM', 'FEKO'}, 'FontSize', 12, 'FontWeight', 'bold');
