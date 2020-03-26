A = load('Spp_conv_N.mat');
Spp = A.Spp;

N = 1:1:30;

figure;

plot(N, db(abs((Spp))), 'LineWidth', 2); grid on;


xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S_{pp} magnitude of TE_{11}'], 'FontSize', 12, 'FontWeight', 'bold')

xlim([1 30]);

