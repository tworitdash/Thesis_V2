clear;

%% Defining radii of each waveguide section
    


rr = 2e-2; % Base radius
rt = 4e-2; % Top redius
n = 15; % number of transitions

R = linspace(rr, rt, n); % radius vector


F = 4e9:0.5e9:21e9; % Frequency of operation

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG section

Length = 5e-2; % height of the cone

L = ones(1, n) .* Length/n; % length of each waveguide section

for i = 1:n
    f =  fc(R(i), er(i), mur(i));
    N_i  =  find(f < F(end));
    N(i) = length(N_i);
end

% N = round(linspace(24, 82, n));

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
STR(k, :, :) = slt * STR_ * slr; 
SRT(k, :, :) = slr * SRT_ * slt; 
SRR(k, :, :) = slr * SRR_ * slr;

end


save('Stt_cone_replica_15_conv', 'STT');
save('Str_cone_replica_15_conv', 'STR');
save('Srt_cone_replica_15_conv', 'SRT');
save('Srr_cone_replica_15_conv', 'SRR');



figure;

% plot(F * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(STT(:, 1, 1))))/2, 'LineWidth', 2); grid on;


figure;

% plot(F * 1e-9, (angle(squeeze(s_params_5(1, 1, :)))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, (angle(squeeze(STT(:, 1, 1)))) * 180/pi, 'LineWidth', 2); grid on;

% %% Plot from files
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
% 
