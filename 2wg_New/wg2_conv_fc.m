A = load('Spp_conv_N_ratio_2_fc.mat');
Spp = A.Spp;

N = 1:1:30;

figure;

plot(N, db(abs((Spp))), 'LineWidth', 2); grid on;


xlabel('Number of modes (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter'], 'FontSize', 12, 'FontWeight', 'bold')

xlim([1 30]);

