
F = linspace(6e9, 12e9, 35);
F1 = linspace(4e9, 21e9, 35);

%% standard mesh


standard = read(rfdata.data,'../../../feko/2wg_lam12.s20p');
S_std = extract(standard,'S_PARAMETERS');


%% Fine mesh 

fine = read(rfdata.data,'../../../feko/2wg_lam20.s20p');
S_fine = extract(fine,'S_PARAMETERS');


%% Super Fine mesh

super_fine = read(rfdata.data,'../../../feko/2wg_lam25.s20p');
S_super_fine = extract(super_fine,'S_PARAMETERS');

%% From MM

c_pp = load('Spp2_ratio_1_modes_20.mat');
SPP = c_pp.Spp;
c_pr = load('Spr2_ratio_1_modes_20.mat');
SPR = c_pr.Spr;
c_rp = load('Srp2_ratio_1_modes_20.mat');
SRP = c_rp.Srp;
c_rr = load('Srr2_ratio_1_modes_20.mat');
SRR = c_rr.Srr;

%% Plots

figure;

plot(F * 1e-9, db(abs(squeeze(S_std(1, 2, :)))), 'LineWidth', 1); grid on;
hold on;

plot(F * 1e-9, db(abs(squeeze(S_fine(1, 2, :)))), 'LineWidth', 1); grid on;
hold on;

plot(F * 1e-9, db(abs(squeeze(S_super_fine(1, 2, :)))), 'LineWidth', 1); grid on;
hold on;


plot(F1 * 1e-9, db(abs(squeeze(SRR(:, 1, 2)))), 'LineWidth', 1); grid on;
hold on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'FEKO \lambda/12', 'FEKO \lambda/20',...
    'FEKO \lambda/25', 'MM'},...
    'FontSize', 12, 'FontWeight', 'bold');

% xlim([9.46 9.48]);
% ylim([-21 -20.5]);
xlim([6 10.5])