clear;

F = 4e9:0.5e9:21e9;


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability


rt = 0.0405319403216/1.9;
rp = 0.0405319403216/2;   % radius of the bigger waveguide
rr = 0.0405319403216/2.1; % radius of the smaller waveguide

R = [rr rp rt];

n = length(R);

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

epsilon = er .* er0;
mu = mur .* mu0;

L = 1e-3 * [1 20 1]; % length of each waveguide sectio

N_ = 30;

for i = 1:5:N_
    
N = [i i i];
    
for k = 1:length(F)

[Spp_, Spr_, Srp_, Srr_] = GSM_N_Conv(R, L, er, mur, F(k), N);

    Spp(i).Spp_i(k, :, :) =  Spp_;
    Spr(i).Spr_i(k, :, :) =  Spr_;
    Srp(i).Srp_i(k, :, :) =  Srp_;
    Srr(i).Srr_i(k, :, :) =  Srr_;

end

% data5 = read(rfdata.data,'S_Feko_5modes_each.s10p');
% s_params_5 = extract(data5,'S_PARAMETERS');

figure(1);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Spp(i).Spp_i(:, 1, 1)))), 'LineWidth', 2); grid on;

figure(2);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Spr(i).Spr_i(:, 1, 1)))), 'LineWidth', 2); grid on;

figure(3);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(1, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Srp(i).Srp_i(:, 1, 1)))), 'LineWidth', 2); grid on;

figure(4);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Srr(i).Srr_i(:, 1, 1)))), 'LineWidth', 2); grid on;

end

save('STT_conv_3', 'Spp');
save('STR_conv_3', 'Spr');
save('SRT_conv_3', 'Srp');
save('SRR_conv_3', 'Srr');

