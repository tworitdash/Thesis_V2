% MM = load('fmincon_SRR.mat');
% SRR = MM.SRR;

MM = load('Paper_minxp_GSM.mat');
SRR_ = MM.b.SRR;
STT_ = MM.b.STT;
STR_ = MM.b.STR;
SRT_ = MM.b.SRT;

MM_ = load('Paper_minxp_GSM_prev.mat');
SRR1 = MM_.b.SRR;
STT1 = MM_.b.STT;
STR1 = MM_.b.STR;
SRT1 = MM_.b.SRT;

% MM_RL = load('fmincon_RL.mat');
% RL = MM_RL.RL;

% F_feko = linspace(4e9, 8e9, 36);
F_feko = linspace(2.6129e9, 7.7561e9, 100);
% F_feko = F;

% data = read(rfdata.data,'../../Thesis/3wg/gsm/RL_fmincon.s2p');
% s_params = extract(data,'S_PARAMETERS');
data = read(rfdata.data,'../../Thesis/3wg/gsm/Paper_Thesis_Closed_minxp_Sparam.s4p');
s_params = extract(data,'S_PARAMETERS');


% RL_feko = db(abs(squeeze(s_params(1, 1, :)) + squeeze(s_params(1, 2, :))).^2)./2;

% figure;
% hold on;
% 
% plot(F * 1e-9, RL, 'LineWidth', 2); hold on;
% plot(F * 1e-9, RL_feko, 'LineWidth', 2); hold on;
% grid on;
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Return loss (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title('RL of ga algo', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'MM', 'FEKO'}, 'FontSize', 12, 'FontWeight', 'bold');

figure;

plot(F_feko * 1e-9, db(abs(squeeze(SRT_(:, 1, 8)))), 'LineWidth', 2); grid on;

hold on;

plot(F * 1e-9, db(abs(squeeze(SRR1(:, 1, 1)))), 'LineWidth', 2); grid on;
hold on;

plot(F * 1e-9, db(abs(squeeze(s_params(3, 3, :)))), 'LineWidth', 2); grid on;
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S11 of TE_{11} with TE_{12} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('S11 in dB', 'FontSize', 12, 'FontWeight', 'bold');
legend({'MM new', 'MM old', 'FEKO'}, 'FontSize', 12, 'FontWeight', 'bold');