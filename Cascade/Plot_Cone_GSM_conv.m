%% 
F_cst = linspace(4e9, 21e9, 1001);
F1 = 1e9:0.5e9:41e9;
F_feko = 5e9:0.5e9:21e9;
F = 4e9:0.5e9:21e9;

% data5_cst =
% read(rfdata.data,'../../../feko/cst/3wg_V2_more_modes.s300p'); % For my mac
% data5_cst = read(rfdata.data,'../../../feko/cst/5wg_V1.s140p'); % For PC at TU
% s_params_5_cst = extract(data5_cst,'S_PARAMETERS');


data5_feko = read(rfdata.data,'../../../feko/Cone_V2.s106p ');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

%%

c_5 = load('Stt_cone_replica_5_conv.mat');
s_5 = c_5.STT;

c_10 = load('Stt_cone_replica_10_conv.mat');
s_10 = c_10.STT;

c_15 = load('Stt_cone_replica_15_conv.mat');
s_15 = c_15.STT;

c_20 = load('Stt_cone_replica_20.mat');
s_20 = c_20.STT;


c_25 = load('Stt_cone_replica_25_conv.mat');
s_25 = c_25.STT;

c_30 = load('Stt_cone_replica_30_conv.mat');
s_30 = c_30.STT;

m_MM = 5;
n_MM = 5;
m = 5;
n = 5;
m_feko = 5;
n_feko = 5;
mode = 'TM';

%% SPP 
figure;
hold on;

plot(F * 1e-9, db(abs(squeeze(s_5(:, m_MM, n_MM)))), 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_10(:, m_MM, n_MM)))), 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_15(:, m_MM, n_MM)))), 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, db(abs(squeeze(s_20(:, m_MM, n_MM)))), 'LineWidth', 2); grid on;
hold on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_25(:, m_MM, n_MM)))), 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_30(:, m_MM, n_MM)))), 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, db(abs(squeeze(s_params_5_feko(m_feko, n_feko, :)))), 'LineWidth', 2); grid on;


xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold')

legend({[mode,'_{', num2str(m),num2str(n),'} MM 5 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 10 elements'],...
    [mode,'_{', num2str(m),num2str(n),'} MM 15 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 20 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 25 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 30 elements'], [mode,'_{', num2str(m),num2str(n),'} FEKO']},...
    'FontSize', 12, 'FontWeight', 'bold');

% xlim([5 9]);
xlim([5 21]);

figure;
hold on;
plot(F * 1e-9, angle((squeeze(s_5(:, m_MM, n_MM)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_10(:, m_MM, n_MM)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_15(:, m_MM, n_MM)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, angle((squeeze(s_20(:, m_MM, n_MM)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_25(:, m_MM, n_MM)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(s_30(:, m_MM, n_MM)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, angle((squeeze(s_params_5_feko(m_feko, n_feko, :)))) * 180/pi, 'LineWidth', 2); grid on;



xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S parameter Phase in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold')

legend({[mode,'_{', num2str(m),num2str(n),'} MM 5 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 10 elements'],...
    [mode,'_{', num2str(m),num2str(n),'} MM 15 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 20 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 25 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 30 elements'], [mode,'_{', num2str(m),num2str(n),'} FEKO']},...
    'FontSize', 12, 'FontWeight', 'bold');



% xlim([5 9]);
xlim([5 21]);

S_5 = abs(squeeze(s_5(:, m_MM, n_MM)));
S_10 = abs(squeeze(s_10(:, m_MM, n_MM)));
S_15 = abs(squeeze(s_15(:, m_MM, n_MM)));
S_20 = abs(squeeze(s_20(:, m_MM, n_MM)));
S_25 = abs(squeeze(s_25(:, m_MM, n_MM)));
S_30 = abs(squeeze(s_30(:, m_MM, n_MM)));
S_feko = abs(squeeze(s_params_5_feko(m_feko, n_feko, :)));
figure;

plot(F_feko * 1e-9, abs(S_feko - S_5(3:end))./S_feko * 100, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, abs(S_feko - S_10(3:end))./S_feko * 100, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, abs(S_feko - S_15(3:end))./S_feko * 100, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, abs(S_feko - S_20(9:41))./S_feko * 100, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, abs(S_feko - S_25(3:end))./S_feko * 100, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, abs(S_feko - S_30(3:end))./S_feko * 100, 'LineWidth', 2); grid on;


legend({[mode,'_{', num2str(m),num2str(n),'} MM 5 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 10 elements'],...
    [mode,'_{', num2str(m),num2str(n),'} MM 15 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 20 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 25 elements'], [mode,'_{', num2str(m),num2str(n),'} MM 30 elements']},...
    'FontSize', 12, 'FontWeight', 'bold');
