spp = load('SPP_conv.mat');

F = 4e9:0.5e9:21e9; % Frequency axis 
figure;


for i = 1:5:30
    
plot(F* 10^(-9), db(abs(spp.Spp(i).Spp_i(:, 1, 1))));
hold on;

end

grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter magnitude in dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TT} of TE_{11} 1 mode', 'S_{TT} of TE_{11} 5 modes',...
    'S_{TT} of TE_{11} 10 modes', 'S_{TT} of TE_{11} 15 modes',...
    'S_{TT} of TE_{11} 20 modes', 'S_{TT} of TE_{11} 25 modes', 'S_{TT} of TE_{11} 30 modes'},...
    'FontSize', 12, 'FontWeight', 'bold');
xlim([5 21]);
