%% 
F_cst = linspace(4e9, 21e9, 1001);
F = 4e9:0.5e9:21e9;
F_feko = 4e9:0.5e9:35e9;

% data5_cst =
% read(rfdata.data,'../../../feko/cst/3wg_V2_more_modes.s300p'); % For my mac
data5_cst = read(rfdata.data,'../../../feko/cst/5wg_V1.s140p'); % For PC at TU
s_params_5_cst = extract(data5_cst,'S_PARAMETERS');


data5_feko = read(rfdata.data,'../../../feko/5wg_V1_feko.s20p');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

%%

% c_pp = load('Stt5_ratio_1_modes_20.mat');
% SPP = c_pp.STT;
% c_pr = load('Str5_ratio_1_modes_20.mat');
% SPR = c_pr.STR;
% c_rp = load('Srt5_ratio_1_modes_20.mat');
% SRP = c_rp.SRT;
% c_rr = load('Srr5_ratio_1_modes_20.mat');
% SRR = c_rr.SRR;


% SPR(:, 3:5, 3:5) = -real(SPR(:, 3:5, 3:5)) - 1j .* imag(SPR(:, 3:5, 3:5));
% SRT(1, 3:5, 3:5) = -real(SRT(1, 3:5, 3:5)) + imag(STR(1, 3:5, 3:5));

%% SPP 
figure;

plot(F * 1e-9, db(abs(squeeze(STT(:, 6, 6)))), 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(9, 9, :)))), 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, db(abs(squeeze(s_params_5_feko(6, 6, :)))), 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TT} of TE_{11} MM', 'S_{TT} of TE_{11} CST',...
    'S_{TT} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

figure;

plot(F * 1e-9, angle((squeeze(STT(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(1, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, angle((squeeze(s_params_5_feko(1, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TT} of TE_{11} MM', 'S_{TT} of TE_{11} CST',...
    'S_{TT} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

%% SPR 

figure;

plot(F * 1e-9, db(abs(squeeze(STR(:, 6, 6)))), 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(8,77, :)))), 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, db(abs(squeeze(s_params_5_feko(6, 16, :)))), 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TR} of TE_{11} MM', 'S_{TR} of TE_{11} CST',...
    'S_{TR} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

figure;

plot(F * 1e-9, angle((squeeze(STR(:, 3, 3)))) * 180/pi, '*', 'LineWidth', 1); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(4, 74, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, angle((squeeze(s_params_5_feko(3, 13, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TR} of TE_{21} MM', 'S_{TR} of TE_{21} CST',...
    'S_{TR} of TE_{21} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

%% SRP

figure;

plot(F * 1e-9, db(abs(squeeze(SRT(:, 5, 5)))), 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(77,7, :)))), 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, db(abs(squeeze(s_params_5_feko(15, 5, :)))), 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RT} of TE_{11} MM', 'S_{RT} of TE_{11} CST',...
    'S_{RT} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);


figure;

plot(F * 1e-9, angle((squeeze(SRT(:, 5, 5)))) * 180/pi, '*', 'LineWidth', 1); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(77, 7, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, angle((squeeze(s_params_5_feko(15, 5, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RT} of TE_{11} MM', 'S_{RT} of TE_{11} CST',...
    'S_{RT} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

%% SRR

figure;

plot(F * 1e-9, db(abs(squeeze(SRR(:, 1, 1)))), 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(71, 71, :)))), 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, db(abs(squeeze(s_params_5_feko(11, 11, :)))), 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RR} of TE_{11} MM', 'S_{RR} of TE_{11} CST',...
    'S_{RR} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);


figure;

plot(F * 1e-9, angle((squeeze(SRR(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(71, 71, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, angle((squeeze(s_params_5_feko(11, 11, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RR} of TE_{11} MM', 'S_{RR} of TE_{11} CST',...
    'S_{RR} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);