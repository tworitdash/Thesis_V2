%% GSM of minxp optimized antenna (matched case)
clear;
% fc_ = fc(R(2), 1, 1);

% F = linspace(fc_(1)+fc_(1)./100, fc_(17), 100);

F = linspace(2.6129e9, 7.7561e9, 100);

MM = load('Paper_minxp_GSM_analyt.mat');
SRR_ = MM.b.SRR;
STT_ = MM.b.STT;
STR_ = MM.b.STR;
SRT_ = MM.b.SRT;

f_int = F(1);

STT = STT_(1, :, :); 
STR = STR_(1, :, :); 

%% Radii of the minxp optimized antenna

data = load('ga_allvar_minxp_V2L.mat');
fmin = data.fmin2;

R = [fmin.r(1) fmin.r(1:5)];

R(end) = 3*R(1);

fc_end = fc(R(end), 1, 1);

[~, indxf] = find(fc_end < f_int);

fc_ = fc(R(1), 1, 1);

[~, indxi] = find(fc_ < f_int);

r = R(end);

%% 
m = 1; % m mode index of the excited mode at the beginning of the conical waveguide

XMN_data = load('Xmn.mat');

xmn = XMN_data.Xmn;
indx_m = 1;
for i = 2:indxf(end)
    if xmn(i).m == m
        indx_m = [indx_m i];
    end
end
indx_mi = 1;
for i = 2:indxi(end)
    if xmn(i).m == m
        indx_mi = [indx_mi i];
    end
end

%% Reform GSM with only the modes having the same azimuthal variation

STT_m = squeeze(STT(1, indx_m, indx_m));
STR_m = squeeze(STR(1, indx_m, indx_mi));

%% Load data of aperture reflection

TE11 = load('Gamma_OpenWG_TE11_freq.mat');
TM11 = load('Gamma_OpenWG_TM11_freq.mat');
TE12 = load('Gamma_OpenWG_TE12_freq.mat');

GTE11 = TE11.Gamma_data.Gamma_11;
GTM11 = TM11.Gamma_data.Gamma_11;
GTE12 = TE12.Gamma_data.Gamma_11;

ar = [1 0];
at = [GTE11 0 0 GTM11 GTE12 0];

% bt = STT_m * at.' + STR_m * ar.';

bt_matched = STR_m * ar.';

%% Find far fields of the modes excited at the final waveguide 

% TE11

dth = pi/180; dphi = pi/180;

[theta, phi] = meshgrid(-pi/2-eps:dth:pi/2, eps:dphi:2*pi+eps);

[Eth_TE11, Eph_TE11, Eco_TE11, Exp_TE11, CO_TE11, XP_TE11] = ...
    FF_TE([1], [1], TE11.Gamma_data.Gamma(1), theta, phi, f_int, 1, 1, r);


[Eth_TM11, Eph_TM11, Eco_TM11, Exp_TM11, CO_TM11, XP_TM11] = ...
    FF_TM([1], [1], TM11.Gamma_data.Gamma_11(1), theta, phi, f_int, 1, 1, r);


[Eth_TE12, Eph_TE12, Eco_TE12, Exp_TE12, CO_TE12, XP_TE12] = ...
    FF_TE([2], [1], TE12.Gamma_data.Gamma(1), theta, phi, f_int, 1, 1, r);

%% For matched waveguide case

[theta, phi] = meshgrid(-pi/2-eps:pi/180:pi/2, eps:pi/180:2*pi+eps);

[Eth_TE11_m, Eph_TE11_m, Eco_TE11_m, Exp_TE11_m, CO_TE11_m, XP_TE11_m] = ...
    FF_TE([1], [1], 0, theta, phi, f_int, 1, 1, r);


[Eth_TM11_m, Eph_TM11_m, Eco_TM11_m, Exp_TM11_m, CO_TM11_m, XP_TM11_m] = ...
    FF_TM([1], [1], 0, theta, phi, f_int, 1, 1, r);


[Eth_TE12_m, Eph_TE12_m, Eco_TE12_m, Exp_TE12_m, CO_TE12_m, XP_TE12_m] = ...
    FF_TE([2], [1], 0, theta, phi, f_int, 1, 1, r);
%% Superposition

% 
Eth = bt_matched(1) .* Eth_TE11 + bt_matched(5) .* Eth_TE12 + bt_matched(4) .* Eth_TM11;
Eph = bt_matched(1) .* Eph_TE11 + bt_matched(5) .* Eph_TE12 + bt_matched(4) .* Eph_TM11;


% Eth = bt_matched(1) .* Eth_TE11;
% Eph = bt_matched(1) .* Eph_TE11;


Eco = bt_matched(1) .* Eco_TE11 + bt_matched(5) .* Eco_TE12 + bt_matched(4) .* Eco_TM11;
Exp = bt_matched(1) .* Exp_TE11 + bt_matched(5) .* Exp_TE12 + bt_matched(4) .* Exp_TM11;


E = sqrt(abs(Eth).^2 + abs(Eph).^2);

Eth_m = bt_matched(1) .* Eth_TE11_m + bt_matched(5) .* Eth_TE12_m + bt_matched(4) .* Eth_TM11_m ;
Eph_m = bt_matched(1) .* Eph_TE11_m + bt_matched(5) .* Eph_TE12_m + bt_matched(4) .* Eph_TM11_m;
% 
% Eth_m = bt_matched(5) .* Eth_TE12_m;
% Eph_m = bt_matched(5) .* Eph_TE12_m;

Eco_m = bt_matched(1) .* Eco_TE11_m + bt_matched(5) .* Eco_TE12_m + bt_matched(4) .* Eco_TM11_m;
Exp_m = bt_matched(1) .* Exp_TE11_m + bt_matched(5) .* Exp_TE12_m + bt_matched(4) .* Exp_TM11_m;


E_m = sqrt(abs(Eth_m).^2 + abs(Eph_m).^2);

% figure; plot(theta(1, :).* 180/pi, db(abs(E_m(1, :))./max(max(abs(E_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(E(1, :))./max(max(abs(E_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% 
% 
% figure; plot(theta(1, :).* 180/pi, db(abs(Eco_m(1, :))./max(max(abs(Eco_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Exp_m(45, :))./max(max(abs(Eco_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Eco(1, :))./max(max(abs(Eco_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Exp(45, :))./max(max(abs(Eco_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);

%% FEKO models
[rho, ph] = meshgrid(linspace(eps, r, 100), linspace(eps, 2*pi, 360));

drho = r/100; dphi = pi/180;

nomfic = '../../Thesis/3wg/gsm/Paper_Minxp_config_NF_OpenWG_2GHz_not_matched_new_cyl.efe';
% nomfic = 'efe';

[E_rho_reshape_, E_phi_reshape_, x_f, y_f] = FEKO_E(nomfic);

[Eth_f, Eph_f, Eco_f, Exp_f] = FF_M_V2(E_rho_reshape_.', E_phi_reshape_.', rho, ph, theta, phi, f_int, drho, dphi);

Ef = sqrt(abs(Eth_f).^2 + abs(Eph_f).^2);

nomfic_m = '../../Thesis/3wg/gsm/Paper_Minxp_config_NF_OpenWG_2GHz_matched_new_cyl.efe';

[E_rho_reshape_m, E_phi_reshape_m, x_fm, y_fm] = FEKO_E(nomfic_m);

[Eth_fm, Eph_fm, Eco_fm, Exp_fm] = FF_M_V2(E_rho_reshape_m.', E_phi_reshape_m.', rho, ph, theta, phi, f_int, drho, dphi);

Ef_m = sqrt(abs(Eth_fm).^2 + abs(Eph_fm).^2);

% figure; plot(theta(1, :).* 180/pi, db(abs(Ef_m(1, :))./max(max(abs(Ef_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Ef(1, :))./max(max(abs(Ef_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% 
% 
% figure; plot(theta(1, :).* 180/pi, db(abs(Eco_fm(1, :))./max(max(abs(Eco_fm)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Exp_fm(45, :))./max(max(abs(Eco_fm)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Eco_f(1, :))./max(max(abs(Eco_fm)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Exp_f(45, :))./max(max(abs(Eco_fm)))), 'LineWidth', 2);  grid on; ylim([-40 0]);




%% 

figure(7); hold on;
% plot(theta(1, :).* 180/pi, db(abs(E_m(1, :))./max(max(abs(E_m(1, :))))), 'LineWidth', 2);  grid on; ylim([-40 0]);
hold on; plot(theta(1, :).* 180/pi, db(abs(E(1, :))./max(max(abs(E_m(1, :))))), 'LineWidth', 2);  grid on; ylim([-40 0]);

% plot(theta(1, :).* 180/pi, db(abs(E_m(1, :))), 'LineWidth', 2);  grid on; ylim([0 35]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(E(1, :))), 'LineWidth', 2);  grid on; %ylim([-40 0]);

% figure(8); hold on;
% plot(theta(1, :).* 180/pi, db(abs(Eph_m(1, :))./max(max(abs(Eph_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Eph(1, :))./max(max(abs(Eph_m)))), 'LineWidth', 2);  grid on; ylim([-40 0]);

hold on;

% plot(theta(1, :).* 180/pi, db(abs(Ef_m(1, :))./max(max(abs(Ef_m)))), '*', 'LineWidth', 1);  grid on; ylim([-40 0]);
% hold on; plot(theta(1, :).* 180/pi, db(abs(Ef(1, :))./max(max(abs(Ef_m)))), '*', 'LineWidth', 1);  grid on; ylim([-40 0]);

% 
% plot(theta(1, :).* 180/pi, db(abs(Ef_m(1, :))), '*', 'LineWidth', 1);  grid on; ylim([0 35]);