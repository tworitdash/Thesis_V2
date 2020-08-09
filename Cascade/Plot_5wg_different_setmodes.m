

F = linspace(4e9, 21e9, 35);
F1 = linspace(4e9, 35e9, 63);

%% From MM

c_pp = load('Stt5_ratio_1_modes_20.mat');
SPP = c_pp.STT;

c_pp_1 = load('Stt5_ratio_1_modes_variable.mat');
SPP_1 = c_pp_1.STT;


data5_feko = read(rfdata.data,'../../../feko/5wg_V1_feko.s20p');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

% c_pr = load('Spr2_ratio_1_modes_20.mat');
% SPR = c_pr.Spr;
% c_rp = load('Srp2_ratio_1_modes_20.mat');
% SRP = c_rp.Srp;
% c_rr = load('Srr2_ratio_1_modes_20.mat');
% SRR = c_rr.Srr;

%% Plots

figure;

plot(F * 1e-9, db(abs(squeeze(SPP(:, 1, 1)))), 'LineWidth', 1); grid on;
hold on;

plot(F * 1e-9, db(abs(squeeze(SPP_1(:, 1, 1)))), 'LineWidth', 1); grid on;
hold on;

% plot(F * 1e-9, db(abs(squeeze(S_super_fine(1, 2, :)))), 'LineWidth', 1); grid on;
% hold on;


plot(F1 * 1e-9, db(abs(squeeze(s_params_5_feko (1, 1, :)))), 'LineWidth', 1); grid on;
hold on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'MM - Same number of modes', 'MM - Different number of modes',  'FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([4 21]);


figure;

plot(F * 1e-9, angle((squeeze(SPP(:, 1, 1)))), 'LineWidth', 1); grid on;
hold on;

plot(F * 1e-9, angle((squeeze(SPP_1(:, 1, 1)))), 'LineWidth', 1); grid on;
hold on;

% plot(F * 1e-9, db(abs(squeeze(S_super_fine(1, 2, :)))), 'LineWidth', 1); grid on;
% hold on;


plot(F1 * 1e-9, angle((squeeze(s_params_5_feko(1, 1, :)))), 'LineWidth', 1); grid on;
hold on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Phase of S in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['Phase of S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'MM - Same number of modes', 'MM - Different number of modes',  'FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([4 21]);