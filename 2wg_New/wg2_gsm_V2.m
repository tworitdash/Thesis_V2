ntic;

clear;

F = 4e9:0.5e9:21e9;

F1 = 4e9:0.5e9:7e9;

rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;

c0 = 3e8;

err = 1;
murr = 1;
erp = 1;
murp = 1;


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;
epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;


Nr = 1:1:20;
Np = 1:1:20;

X_til = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr);


parfor k = 1:length(F)


[Spp_, Spr_, Srp_, Srr_] = GSM(Nr, Np, F(k), rp, rr, erp, murp, err, murr, X_til);

% Spp(k, :, :) = Spp_;
% Spr(k, :, :) = Spr_;
% Srp(k, :, :) = Srp_;
% Srr(k, :, :) = Srr_;

Slp = SL(rr, F(k), Np, -0.001);
Slr = SL(rp, F(k), Nr, 0.001);

Spp(k, :, :) = Slp * Spp_ * Slp';
Spr(k, :, :) = Slp * Spr_ * Slr';
Srp(k, :, :) = Slr * Srp_ * Slp';
Srr(k, :, :) = Slr * Srr_ * Slr';

S = [Spp_ Spr_; Srp_ Srr_];
S*S


end



% save('Spp2_ratio_1_modes_20_test', 'Spp');
% save('Spr2_ratio_1_modes_20_test', 'Spr');
% save('Srp2_ratio_1_modes_20_test', 'Srp');
% save('Srr2_ratio_1_modes_20_test', 'Srr');


%%
% F_cst = linspace(4e9, 21e9, 1001);
% data5_cst = read(rfdata.data,'2wg_cst.s38p');
% s_params_5_cst = extract(data5_cst,'S_PARAMETERS');
% 
% 
% data5_feko = read(rfdata.data,'2wg_feko.s20p');
% s_params_5_feko = extract(data5_feko,'S_PARAMETERS');
% 
% figure;
% 
% plot(F * 1e-9, db(abs(squeeze(Spr(:, 1, 5))))/2, 'LineWidth', 2); grid on;
% % hold on;
% % plot(F * 1e-9, db(abs(squeeze(Srp(:, 1, 5))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(1,25, :))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(1, 15, :))))/2, 'LineWidth', 2); grid on;
% 
% 
% xlim([5 21])
% 
% % figure;
% 
% figure;
% 
% plot(F * 1e-9, angle((squeeze(Spr(:, 1, 5)))), 'LineWidth', 2); grid on;
% % hold on;
% % plot(F * 1e-9, angle((squeeze(Srp(:, 1, 5)))), 'LineWidth', 2); grid on;
% hold on;
% plot(F_cst * 1e-9, angle((squeeze(s_params_5_cst(1,25, :)))), 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, angle((squeeze(s_params_5_feko(5, 11, :)))), 'LineWidth', 2); grid on;


%% 

% figure;
% 
% hold on;
% 
% plot(F * 1e-9, db(abs(squeeze(Srr(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% hold on;
% 
% figure;
% 
% plot(F * 1e-9, db(abs(squeeze(Spr(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(Srp(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% figure;
% 
% plot(F * 1e-9, db(abs(squeeze(Srr(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% hold on;

wg2_time = toc;