%% 
F_cst = linspace(4e9, 21e9, 1001);
F = 4e9:0.5e9:21e9;

% data5_cst =read(rfdata.data,'3wg_V1.s39p'); % For my mac
data5_cst = read(rfdata.data,'/tudelft.net/staff-groups/ewi/me/MS3/Studenten/Tworit Dash/Feko_Files/feko/cst/3wg_V2_more_modes.s300p'); % For PC at TU
s_params_5_cst = extract(data5_cst,'S_PARAMETERS');


data5_feko = read(rfdata.data,'../Cascade/3wg_feko_10_modes.s20p');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

%%
% 
% c_pp = load('../Cascade/Stt3_ratio_1_modes_20_1mm.mat');
% SPP = c_pp.STT;
% c_pr = load('../Cascade/Str3_ratio_1_modes_20_1mm.mat');
% SPR = c_pr.STR;
% c_rp = load('../Cascade/Srt3_ratio_1_modes_20_1mm.mat');
% SRP = c_rp.SRT;
% c_rr = load('../Cascade/Srr3_ratio_1_modes_20_1mm.mat');
% SRR = c_rr.SRR;

%%
c_pp = load('Stt3_ratio_1_modes_20_new_analyt.mat');
SPP = c_pp.Stt;
c_pr = load('Str3_ratio_1_modes_20_new_analyt.mat');
SPR = c_pr.Str;
c_rp = load('Srt3_ratio_1_modes_20_new_analyt.mat');
SRP = c_rp.Srt;
c_rr = load('Srr3_ratio_1_modes_20_new_analyt.mat');
SRR = c_rr.Srr;


%%

% c_pp = load('Stt3_ratio_1_modes_20_new.mat');
% SPP = c_pp.Stt;
% c_pr = load('Str3_ratio_1_modes_20_new.mat');
% SPR = c_pr.Str;
% c_rp = load('Srt3_ratio_1_modes_20_new.mat');
% SRP = c_rp.Srt;
% c_rr = load('Srr3_ratio_1_modes_20_new.mat');
% SRR = c_rr.Srr;
%% SPP 
figure;

plot(F * 1e-9, db(abs(squeeze(SPP(:, 1, 1)))), 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(1,1, :)))), 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(11, 11, :)))), 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TT} of TE_{11} MM', 'S_{TT} of TE_{11} CST',...
    'S_{TT} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

figure;

plot(F * 1e-9, angle((squeeze(SPP(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(1, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_params_5_feko(11, 11, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TT} of TE_{11} MM', 'S_{TT} of TE_{11} CST',...
    'S_{TT} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

%% SPR 

figure;

plot(F * 1e-9, db(abs(squeeze(SPR(:, 1, 8)))), 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(1,8, :)))), 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(1, 6, :)))), 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TR} of TE_{11} MM', 'S_{TR} of TE_{11} CST',...
    'S_{TR} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

figure;

plot(F * 1e-9, angle((squeeze(SPR(:, 1, 1)))) * 180/pi, '*', 'LineWidth', 1); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(1, 151, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_params_5_feko(1, 11, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TR} of TE_{11} MM', 'S_{TR} of TE_{11} CST',...
    'S_{TR} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

%% SRP

figure;

plot(F * 1e-9, db(abs(squeeze(SRP(:, 1, 1)))), 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(71,1, :)))), 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(6, 1, :)))), 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RT} of TE_{11} MM', 'S_{RT} of TE_{11} CST',...
    'S_{RT} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);


figure;

plot(F * 1e-9, angle((squeeze(SRP(:, 1, 1)))) * 180/pi, '*', 'LineWidth', 1); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(71, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_params_5_feko(6, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;

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
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(151, 151, :)))), 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(1, 1, :)))), 'LineWidth', 2); grid on;

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
plot(F * 1e-9, angle((squeeze(s_params_5_feko(6, 6, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RR} of TE_{11} MM', 'S_{RR} of TE_{11} CST',...
    'S_{RR} of TE_{11} FEKO'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);