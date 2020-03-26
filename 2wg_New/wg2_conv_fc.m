A = load('Spr_conv_N_ratio_2.mat');
Spp = A.Spr;

N = 1:1:30;

figure;

plot(N, db(abs((Spp))), 'LineWidth', 2); grid on;


xlabel('Number of modes (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S_{pp} magnitude of TE_{11}'], 'FontSize', 12, 'FontWeight', 'bold')

xlim([1 30]);

