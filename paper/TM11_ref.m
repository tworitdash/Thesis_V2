clear;
tic;

c0 = 3e8; % speed of light
% R = 2e-2; % Radius of the waveguide

er = 1;
mur = 1;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er .* er0;
mu = mur .* mu0;
% f = linspace(10.15e9, 14e9, 41);
f = linspace(18.597e9, 30e9, 41);
%f = linspace(14.1569e9, 16.6551e9, 41);

% f = linspace(22.7e9, 25e9, 20);
% f = linspace(14.1569e9, 70.6551e9, 200);
% f = linspace(23e9, 30e9, 20);
% f = linspace(4.99654e9, 8.32757e9, 41);
% f = linspace(4.99654e9, 20e9, 41);
n_orig = 2;
% F = 5e9;
% lamb = c0/F;

% x = 0.6:0.01:1;
% x = 1.7:0.01:2;
% x = linspace(1.7, 2, 41);
% r = x .* lamb ./ 2;
r = 1.8e-2;
% F = 5e9;

% lamb = c0/F;
omega = 2 * pi * f;

% len = (2 .* r(1))./4; % To match with FEKO
len = 0.036/4;

k0 = omega./c0;

L = 100;

N = 3:1:6;
% N = [3 4 5];

[Dm_, yap, Gamma, y11, Gamma_11, YTE11, y_for_debug] = Tworit_Integrals_K_Space_freq_TM(r, N, k0, L, er, mur, n_orig);

% Gamma = (1 - yap)./(yap + 1);
% Gamma_11 = (1 - y11)./(y11 + 1);
% Gamma = Gamma'

% figure(4);
% hold on;
% plot(f, real(YTE11));
% hold on;
% plot(f, real(y_for_debug));
% hold on;
% plot(f, ones(length(YTE11)) .* 1./(c0 .* mu), 'g');

time = toc;

% 
% [Yin, Gamma_K] = Mishustin_K_Freq_TM(r, k0, L, er, mur, len, n_orig);

%% FEKO results
% 
% data5_feko = read(rfdata.data,'/home/nfs/tworitdash/tdash/Thesis/Paper_Thesis/WG_1_1_8_GP.s1p');
% s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

% data5_feko_GP10r = read(rfdata.data,'/home/nfs/tworitdash/tdash/Thesis/Paper_Thesis/WG_1_1_8_5GHz_10rGP.s3p');
% data5_feko_GP10r = read(rfdata.data,'/home/nfs/tworitdash/tdash/Thesis/Paper_Thesis/WG_1_1_8_GP_TE12_Big_Big_GP.s1p');
data5_feko_GP10r = read(rfdata.data,'/home/nfs/tworitdash/tdash/Thesis/Paper_Thesis/WG_1_1_8_GP_TM11.s1p');
s_params_5_feko = extract(data5_feko_GP10r,'S_PARAMETERS');




m_feko = 1;
n_feko = 1;

Str = load('Xmn_azimuthal_inc_TM.mat');
% f_ = f;
% f_ = linspace(4.99654e9, 8.32757e9, 41);
% r_feko = 1.8e-2;
% f_ = linspace(14.1569e9, 16.6551e9, 41);
f_ = f;
% f_ = linspace(23e9, 30e9, 20);

k0_ = (2 * pi * f_)/c0;

str = Str.xmn_TM;

% xmn11 = str(n_feko).xmn;
xmn11 = str(1).xmn;
beta_rho11 = xmn11./r;

beta_z11 = -1j .* sqrt(-(k0_.^2 - beta_rho11.^2));

Gamma_FEKO = squeeze(s_params_5_feko(m_feko, n_feko, :)).' .* exp(+2j .* beta_z11 .* len);
Y_FEKO = (1 - Gamma_FEKO)./(1 + Gamma_FEKO);

% 
% Gamma_FEKO_Gp10r = squeeze(s_params_5_feko_GP10r(m_feko, n_feko, :)).'  .* exp(+2j .* beta_z11 .* len);
% Y_FEKO_Gp10r = (1 - Gamma_FEKO_Gp10r)./(1 + Gamma_FEKO_Gp10r);

%% CST Results
% 
% data5_feko = read(rfdata.data,'/home/nfs/tworitdash/tdash/Thesis/Paper_Thesis/WG_1_1_8_GP.s1p');
% s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

% data5_cst_GP10r = read(rfdata.data,'/home/nfs/tworitdash/tdash/Thesis/Paper_Thesis/WG_ref_FS_TE11.s1p');

data5_cst_GP10r = read(rfdata.data,'/home/nfs/tworitdash/tdash/Thesis/Paper_Thesis/WG_ref_FS_TE12.s20p');
s_params_5_cst = extract(data5_cst_GP10r,'S_PARAMETERS');


m_cst = 16;
n_cst = 16;

Str = load('Xmn_azimuthal_inc_TE.mat');
% f_cst = linspace(5e9, 20e9, 1001);
f_cst = linspace(14.1569e9, 16.6551e9, 1001);
% r_feko = 1.8e-2;
% f_ = linspace(14.1569e9, 16.6551e9, 41);
% f_ = linspace(23e9, 30e9, 20);

k0_cst = (2 * pi * f_cst)/c0;

str = Str.xmn_TE;

% xmn11 = str(n_feko).xmn;
xmn11 = str(2).xmn;
beta_rho11 = xmn11./r;

beta_z11_cst = -1j .* sqrt(-(k0_cst.^2 - beta_rho11.^2));

Gamma_CST = squeeze(s_params_5_cst(m_cst, n_cst, :)).' .* exp(2j .* beta_z11_cst .* len);
Y_CST = (1 - Gamma_CST)./(1 + Gamma_CST);
%%

figure;
hold on;
plot(f*1e-9, real(yap), 'o', 'Linewidth', 1);
hold on;
plot(f*1e-9, imag(yap), 'o', 'Linewidth', 1);

grid on;
hold on;
plot(f*1e-9, real(y11), '-.', 'Linewidth', 1);
hold on;
plot(f*1e-9, imag(y11), '-.', 'Linewidth', 1);

% hold on;
% plot(f*1e-9, real(Yin), '*', 'Linewidth', 1);
% hold on;
% plot(f*1e-9, imag(Yin), '*', 'Linewidth', 1);


% 
% hold on;
% plot(x, real(Y_FEKO_Gp10r), '+', 'Linewidth', 1);
% hold on;
% plot(x, imag(Y_FEKO_Gp10r), '+', 'Linewidth', 1);
hold on;
plot(f_*1e-9, real(Y_FEKO), 'Linewidth', 1);
hold on;
plot(f_*1e-9, imag(Y_FEKO), 'Linewidth', 1);
% hold on;
% plot(f_cst*1e-9, real(Y_CST), 'Linewidth', 1);
% hold on;
% plot(f_cst*1e-9, imag(Y_CST), 'Linewidth', 1);

% xlim([4.99654 8.32757]);
% xlim([14.1569 16.6551])
% xlim([22.7, 23.426])
xlabel('Frequency (GHz)', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('y_{ap}', 'FontWeight', 'bold', 'FontSize', 16);
title('Aperture Admittance over Free space admittance', 'FontWeight', 'bold', 'FontSize', 16);

% legend({'Real y_{ap} This technique with 3 higher order modes', ...
%     'Imag y_{ap} This technique with 3 higher order modes', ...
%     'Real y_{ap} This technique with no higher order modes'...
%     ,'Imag y_{ap} This technique with no higher order modes', ...
%     'Real y_{ap} Mishustin'...
%     ,'Imag y_{ap} Mishustin', ...
%     'Real y_{ap} FEKO', 'Imag y_{ap} FEKO',...
%     'Real y_{ap} CST', 'Imag y_{ap} CST'})

legend({'Real y_{ap} This technique with 3 higher order modes', ...
    'Imag y_{ap} This technique with 3 higher order modes', ...
    'Real y_{ap} This technique with no higher order modes'...
    ,'Imag y_{ap} This technique with no higher order modes', ...
    'Real y_{ap} Mishustin'...
    ,'Imag y_{ap} Mishustin', ...
    'Real y_{ap} FEKO', 'Imag y_{ap} FEKO'})

%%
figure;
hold on;
plot(f*1e-9, db(abs(Gamma)), 'o', 'Linewidth', 1);
hold on;
plot(f*1e-9, db(abs(Gamma_11)), '-.', 'Linewidth', 1);
hold on;
% plot(f*1e-9, db(abs(Gamma_K)), '*', 'Linewidth', 1);
grid on;
hold on;
plot(f_*1e-9, db(abs(Gamma_FEKO)), 'r', 'LineWidth', 2); grid on;
hold on;
% plot(f_cst*1e-9, db(abs(Gamma_CST)), 'k', 'LineWidth', 2); grid on;
% hold on;
% plot(x, db(abs(Gamma_FEKO_Gp10r)), '+', 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('\Gamma in dB', 'FontWeight', 'bold', 'FontSize', 16);
title('Reflection Coefficient', 'FontWeight', 'bold', 'FontSize', 16);
% xlim([4.99654 8.32757]);
% xlim([14.1569 16.6551])
% xlim([22.7, 23.426])
% 
% legend({'\Gamma This technique with 3 higher order modes', ...
%     '\Gamma This technique with no higher order modes', ...
%     '\Gamma Mishustin',...
%     '\Gamma FEKO',...
%     '\Gamma CST'});

legend({'\Gamma This technique with 3 higher order modes', ...
    '\Gamma This technique with no higher order modes', ...
    '\Gamma Mishustin',...
    '\Gamma FEKO'});
%%
figure;
hold on;
plot(f, unwrap(angle((Gamma))) .* 180/pi, 'o', 'Linewidth', 1);
hold on;
plot(f, angle((Gamma_11)) .* 180/pi, '-.', 'Linewidth', 1);
hold on;
% plot(f, angle((Gamma_K)) .* 180/pi, '*', 'Linewidth', 1);
grid on;
hold on;
plot(f_, unwrap(angle((Gamma_FEKO))) .* 180/pi, 'r', 'LineWidth', 2); grid on;

hold on;    
% plot(f_cst, -unwrap(angle((Gamma_CST))) .* 180/pi, 'k', 'LineWidth', 2); grid on;
% hold on;
% plot(x, db(abs(Gamma_FEKO_Gp10r)), '+', 'LineWidth', 2); grid on;

xlabel('2r/ \lambda', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Phase of \Gamma in deg', 'FontWeight', 'bold', 'FontSize', 16);
title('Reflection Coefficient', 'FontWeight', 'bold', 'FontSize', 16);


