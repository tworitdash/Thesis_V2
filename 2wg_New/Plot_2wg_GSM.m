
%% 
F_cst = linspace(4e9, 21e9, 1001);
F = 4e9:0.5e9:21e9;

data5_cst = read(rfdata.data,'2wg_cst.s38p');
s_params_5_cst = extract(data5_cst,'S_PARAMETERS');


data5_feko = read(rfdata.data,'2wg_feko.s20p');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

c_pp = load('Spp2_ratio_1_modes_20.mat');
SPP = c_pp.Spp;
c_pr = load('Spr2_ratio_1_modes_20.mat');
SPR = c_pr.Spr;
c_rp = load('Srp2_ratio_1_modes_20.mat');
SRP = c_rp.Srp;
c_rr = load('Srr2_ratio_1_modes_20.mat');
SRR = c_rr.Srr;

%% SPP 
figure;

plot(F * 1e-9, db(abs(squeeze(SPP(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(1,1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(11, 11, :))))/2, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{PP} of TE_{11} 20 modes active MM', 'S_{PP} of TE_{11} 10 modes active CST',...
    'S_{PP} of TE_{11}, 10 modes active MM'},...
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

legend({'S_{PP} of TE_{11} 20 modes active MM', 'S_{PP} of TE_{11} 10 modes active CST',...
    'S_{PP} of TE_{11}, 10 modes active MM'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

%% SPR 

figure;

plot(F * 1e-9, db(abs(squeeze(SPR(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(1,20, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(11, 1, :))))/2, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{PR} of TE_{11} 20 modes active MM', 'S_{PR} of TE_{11} 10 modes active CST',...
    'S_{PR} of TE_{11}, 10 modes active MM'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

figure;

plot(F * 1e-9, angle((squeeze(SPR(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(1, 20, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_params_5_feko(11, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{PR} of TE_{11} 20 modes active MM', 'S_{PR} of TE_{11} 10 modes active CST',...
    'S_{PR} of TE_{11}, 10 modes active MM'},...
    'FontSize', 12, 'FontWeight', 'bold');


xlim([5 21]);

%% SRP

figure;

plot(F * 1e-9, db(abs(squeeze(SRP(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(20,1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(1, 11, :))))/2, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RP} of TE_{11} 20 modes active MM', 'S_{RP} of TE_{11} 10 modes active CST',...
    'S_{RP} of TE_{11}, 10 modes active MM'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);


figure;

plot(F * 1e-9, angle((squeeze(SRP(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(20, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_params_5_feko(1, 11, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RP} of TE_{11} 20 modes active MM', 'S_{RP} of TE_{11} 10 modes active CST',...
    'S_{RP} of TE_{11}, 10 modes active MM'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

%% SRR

figure;

plot(F * 1e-9, db(abs(squeeze(SRR(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(20,20, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(1, 1, :))))/2, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RR} of TE_{11} 20 modes active MM', 'S_{RR} of TE_{11} 10 modes active CST',...
    'S_{RR} of TE_{11}, 10 modes active MM'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);


figure;

plot(F * 1e-9, angle((squeeze(SRR(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(20, 20, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_params_5_feko(1, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{RR} of TE_{11} 20 modes active MM', 'S_{RR} of TE_{11} 10 modes active CST',...
    'S_{RR} of TE_{11}, 10 modes active MM'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([5 21]);

