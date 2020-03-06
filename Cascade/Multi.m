clear;

%% Defining radii of each waveguide section
    
rt = 0.0405319403216/1.9;
rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;
rd = 2.2e-2;
re = 2.3e-2;


R = [rr rp rt rd re]; % radius vector

F = 4e9:0.5e9:21e9; % Frequency of operation

er = [1 1 1 1 1]; % Relative Permittivity of each WG section
mur = [1 1 1 1 1]; % Relative Permeability of each WG section

L = 1e-3 * [1 1 1 20 1]; % length of each waveguide section

N = 1:1:20; % Number of modes

J = length(R) - 1; % Number of Junctions
%% Frequency independent inner cross product 

for j = 1:J

    X_til(j, :, :) = Inner_p(N, N, R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j));

end

%% Frequency loop to find the GSM of the entire structure

parfor k = 1:length(F)
    
    disp('Frequency Iteration: ');
    disp(k);

[S33, S34, S43, S44] = GSM(N, N, F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), squeeze(X_til(1, :, :)));
[S11, S12, S21, S22] = GSM(N, N, F(k), R(3), R(2), er(3), mur(3), er(2), mur(2), squeeze(X_til(2, :, :)));
Sl = SL(R(2), F(k), N, L(2));

  
[STT_, STR_, SRT_, SRR_] = cascade_3(N, S11, S12, S21, S22, S33, S34, S43, S44, Sl);

% Use the for loop in case of more than 3 junctions (J > 3)

if J > 2

for j = 3:J

    % recursion 
    
    [S11, S12, S21, S22] = GSM(N, N, F(k), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j), squeeze(X_til(j, :, :)));
    S33 = STT_; S34 = STR_; S43 = SRT_; S44 = SRR_;
    Sl = SL(R(j), F(k), N, L(j));
    
    [STT_, STR_, SRT_, SRR_] = cascade_3(N, S11, S12, S21, S22, S33, S34, S43, S44, Sl);
    
end

else 
    
end

slr = SL(R(1), F(k), N, L(1));
slt = SL(R(end), F(k), N, L(end));

STT(k, :, :) = slt * STT_ * slt'; 
STR(k, :, :) = slt * STR_ * slr; 
SRT(k, :, :) = slr * SRT_ * slt; 
SRR(k, :, :) = slr * SRR_ * slr;

end


save('Stt5_ratio_1_modes_20', 'STT');
save('Str5_ratio_1_modes_20', 'STR');
save('Srt5_ratio_1_modes_20', 'SRT');
save('Srr5_ratio_1_modes_20', 'SRR');



%% Plot


data5 = read(rfdata.data,'5wg_touchstone_5modes_1mm_2cm_1mm.s10p');
s_params_5 = extract(data5,'S_PARAMETERS');

c_ = load('Stt5_ratio_1_modes_20.mat');

F1 = 4e9:0.5e9:21e9;

STT = c_.STT;
figure;
plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(STT(:, 1, 1))))/2, 'LineWidth', 2); grid on;
xlim([4.5 21]);
figure;
plot(F1 * 1e-9, angle((squeeze(s_params_5(1, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, (angle(squeeze(STT(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;

xlim([4.5 21]);

c_ = load('Str5_ratio_1_modes_20.mat');

STR = c_.STR;
figure;
plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(STR(:, 1, 1))))/2, 'LineWidth', 2); grid on;
xlim([4.5 21]);
figure;
plot(F1 * 1e-9, angle((squeeze(s_params_5(1, 6, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, (angle(squeeze(STR(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;

xlim([4.5 21]);

c_ = load('Srt5_ratio_1_modes_20.mat');
SRT = c_.SRT;
figure;
plot(F1 * 1e-9, db(abs(squeeze(s_params_5(6, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(SRT(:, 1, 1))))/2, 'LineWidth', 2); grid on;
xlim([4.5 21]);
figure;
plot(F1 * 1e-9, angle((squeeze(s_params_5(6, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, (angle(squeeze(SRT(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;

xlim([4.5 21]);


c_ = load('Srr5_ratio_1_modes_20.mat');
SRR = c_.SRR;
figure;
plot(F1 * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(SRR(:, 1, 1))))/2, 'LineWidth', 2); grid on;
xlim([4.5 21]);
figure;
plot(F1 * 1e-9, angle((squeeze(s_params_5(6, 6, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, (angle(squeeze(SRR(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;

xlim([4.5 21]);

