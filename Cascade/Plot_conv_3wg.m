spp = load('STT_conv_3.mat');

F = 4e9:0.5e9:21e9; % Frequency axis 
figure;


for i = [1  11 21 26]
    
    if i == 26
        str = '*';
    else
        str = '-';
    end
    
plot(F* 10^(-9), db(abs(spp.Spp(i).Spp_i(:, 1, 1))), str, 'Color', [1, 1 - (i + 10)/40 , 1], 'LineWidth', 2);
hold on;

end

grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter magnitude in dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({'S_{TT} of TE_{11} 1 mode',...
    'S_{TT} of TE_{11} 11 modes',...
    'S_{TT} of TE_{11} 21 modes', 'S_{TT} of TE_{11} 26 modes'},...
    'FontSize', 12, 'FontWeight', 'bold');
xlim([5 21]);
