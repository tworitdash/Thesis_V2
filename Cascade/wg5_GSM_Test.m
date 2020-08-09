clear;

%% Defining radii of each waveguide section
    

rt = 0.0405319403216/1.9;
rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;
rd = 2.2e-2;
re = 2.3e-2;


R = [rr rp rt rd re]; % radius vector

n = length(R);


F = 4e9:0.5e9:21e9; % Frequency of operation
F_feko = F;

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG section

% Length = 5e-2; % height of the cone

% L = ones(1, n) .* Length/n; % length of each waveguide section


L = 1e-3 * [1 1 1 20 1]; % length of each waveguide section

[fc_1] = fc(R(1), er(1), mur(1));
[fc_2] = fc(R(2), er(2), mur(2));
[fc_3] = fc(R(3), er(3), mur(3));
[fc_4] = fc(R(4), er(4), mur(4));
[fc_5] = fc(R(5), er(5), mur(5));

N1 = find(fc_1 < F(end));
N2 = find(fc_2 < F(end));
N3 = find(fc_3 < F(end));
N4 = find(fc_4 < F(end));
N5 = find(fc_5 < F(end));

N = [N1(end) N2(end) N3(end) N4(end) N5(end)];

% N = round(linspace(5, 25, n));

% N = 1:1:5 ; % Number of modes

J = length(R) - 1; % Number of Junctions
%% Frequency independent inner cross product 

for j = 1:J
    x_til = zeros(N(j), N(j + 1));
    x_til(:, :) = Inner_p(1:1:N(j), 1:1:N(j + 1), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j));
    X_til(j).x_til = x_til;
end

%% Frequency loop to find the GSM of the entire structure

parfor k = 1:length(F)
    
    disp('Frequency Iteration: ');
    disp(k);

[S33, S34, S43, S44] = GSM(1:1:N(1), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), X_til(1).x_til);
[S11, S12, S21, S22] = GSM(1:1:N(2), 1:1:N(3), F(k), R(3), R(2), er(3), mur(3), er(2), mur(2), X_til(2).x_til);
Sl = SL(R(2), F(k), 1:1:N(2), L(2));

  
[STT_, STR_, SRT_, SRR_] = cascade_3(1:1:N(2), S11, S12, S21, S22, S33, S34, S43, S44, Sl);

% Use the for loop in case of more than 3 junctions (J > 3)

for j = 3:J

    % recursion 
    
    [S11, S12, S21, S22] = GSM(1:1:N(j), 1:1:N(j + 1), F(k), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j), X_til(j).x_til);
    S33 = STT_; S34 = STR_; S43 = SRT_; S44 = SRR_;
    Sl = SL(R(j), F(k), 1:1:N(j), L(j));
    
    [STT_, STR_, SRT_, SRR_] = cascade_3(1:1:N(j), S11, S12, S21, S22, S33, S34, S43, S44, Sl);
    
end

slr = SL(R(1), F(k), 1:1:N(1), L(1));
slt = SL(R(end), F(k), 1:1:N(end), L(end));

STT(k, :, :) = slt * STT_ * slt'; 
STR(k, :, :) = slt * STR_ * slr'; 
SRT(k, :, :) = slr * SRT_ * slt'; 
SRR(k, :, :) = slr * SRR_ * slr';

end


save('Stt5_ratio_1_modes_variable_test', 'STT');
save('Str5_ratio_1_modes_variable_test', 'STR');
save('Srt5_ratio_1_modes_variable_test', 'SRT');
save('Srr5_ratio_1_modes_variable_test', 'SRR');

data5_feko = read(rfdata.data,'../../../feko/5wg_V1_feko.s20p');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

figure;

% plot(F * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(STR(:, 3, 3))))/2, 'LineWidth', 2); grid on;
hold on;
hold on;
plot(F_feko * 1e-9, db(abs(squeeze(s_params_5_feko(3, 13, :)))), 'LineWidth', 2); grid on;

figure;

% plot(F * 1e-9, (angle(squeeze(s_params_5(1, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, (angle(squeeze(STR(:, 3, 3)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F_feko * 1e-9, angle((squeeze(s_params_5_feko(3, 13, :)))) * 180/pi, 'LineWidth', 2); grid on;


%% Plot from files
% % 
% % data5 = read(rfdata.data,'Cone_feko.s2p');
% % s_params_5 = extract(data5,'S_PARAMETERS');
% 
% c_ = load('Stt_cone_replica_20.mat');
% 
% F1 = 1e9:0.5e9:21e9;
% 
% STT = c_.STT;
% figure;
% % plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(STT(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% % figure(2);
% % % plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
% % hold on;
% % plot(F * 1e-9, (angle(squeeze(STT(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
% 
% xlim([4.5 21]);
% 
% c_ = load('Str_cone_replica_20.mat');
% 
% STR = c_.STR;
% figure;
% % plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(STR(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% % figure(2);
% % % plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
% % hold on;
% % plot(F * 1e-9, (angle(squeeze(STT(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
% 
% xlim([4.5 21]);
% 
% c_ = load('Srt_cone_replica_20.mat');
% SRT = c_.SRT;
% 
% % plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(SRT(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% % figure(2);
% % % plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
% % hold on;
% % plot(F * 1e-9, (angle(squeeze(STT(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
% 
% xlim([4.5 21]);
% 
% 
% c_ = load('Srr_cone_replica_20.mat');
% SRR = c_.SRR;
% figure;
% % plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(SRR(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% % figure(2);
% % % plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
% % hold on;
% % plot(F * 1e-9, (angle(squeeze(STT(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;
% 
% xlim([4.5 21]);

