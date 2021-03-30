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

%f = linspace(14.1569e9, 16.6551e9, 41);

% f = linspace(22.7e9, 25e9, 20);

% f = linspace(14.1569e9, 70.6551e9, 200);
% f = linspace(23e9, 30e9, 20);
% f = linspace(4.99654e9, 8.32757e9, 41);
f = 17e9;
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

N = 3:1:10;
% N = [3 4 5];

[Dm_, yap, Gamma, y11, Gamma_11, YTE11, y_for_debug] = Tworit_Integrals_K_Space_freq(r, N, k0, L, er, mur, n_orig);
% Gamma = Gamma'

% figure(4);
% hold on;
% plot(f, real(YTE11));
% hold on;
% plot(f, real(y_for_debug));
% hold on;
% plot(f, ones(length(YTE11)) .* 1./(c0 .* mu), 'g');

time = toc;

% % 
% [Yin, Gamma_K] = Mishustin_K_Freq_TE12(r, k0, L, er, mur, len, n_orig);

%% Near Fields
[rho, phi] = meshgrid(linspace(eps, r, 100), linspace(eps, 2*pi, 360));
Str = load('Xmn_azimuthal_inc_TE.mat');
Xmn = Str.xmn_TE;
xmn11 = Xmn(2).xmn;
beta_rho11 = xmn11/r;

Dmf = Dm_;

[Erho_f, Ephi_f, Ez_f, Hrho_f, Hphi_f, Hz_f] = E_TE_FEKOstyle(mu, 1, rho, phi, beta_rho11, len, k0(1));
Ef_nr = sqrt(abs(Erho_f).^2 + abs(Ephi_f).^2);
Hf_nr = sqrt(abs(Hrho_f).^2 + abs(Hphi_f).^2);
Erho = (1 + Gamma(1)) * Erho_f;
Ephi = (1 + Gamma(1)) * Ephi_f;

Hrho = (1 - Gamma(1)) * Hrho_f;
Hphi = (1 - Gamma(1)) * Hphi_f;

Ef = sqrt(abs(Erho).^2 + abs(Ephi).^2);
Hf = sqrt(abs(Hrho).^2 + abs(Hphi).^2);

  figure(1);
  plot(rho(1, :), Ef_nr(1, :).'); grid on;
  hold on;
  plot(rho(1, :), Ef(1, :).'); grid on;
   
  figure(2);
  plot(rho(1, :), Hf_nr(1, :).'); grid on;
  hold on;
  plot(rho(1, :), abs(Hf(1, :).')); grid on;
% 
% figure;
% surface(rho .* cos(phi), rho .* sin(phi), abs(Ef)); shading flat; colormap('jet');
% 
% figure;
% surface(rho .* cos(phi), rho .* sin(phi), abs(Hf));  shading flat; colormap('jet');

for k = 2:length(N)
 
  xmn = Xmn(N(k)).xmn;
  beta_rho = xmn/r;
  [Erho_i, Ephi_i, Ez_i, Hrho_i, Hphi_i, Hz_i] = E_TE_FEKOstyle(mu, 1, rho, phi, beta_rho, len, k0(1));
  Erho = Erho + (1 + Gamma(1)) .* Erho_i * Dmf(k);
  Ephi = Ephi + (1 + Gamma(1)) .* Ephi_i * Dmf(k);
  Hrho = Hrho + (-1 - Gamma(1)) .* Hrho_i * Dmf(k);
  Hphi = Hphi + (-1 - Gamma(1)) .* Hphi_i * Dmf(k);
  
  Ef = sqrt(abs(Erho).^2 + abs(Ephi).^2);
  Hf = sqrt(abs(Hrho).^2 + abs(Hphi).^2);

  figure(1);
  hold on;
  plot(rho(1, :), Ef(1, :));
  figure(2);
  hold on;
  plot(rho(1, :), abs(Hf(1, :)));
  
%   figure;
%   surface(rho .* cos(phi), rho .* sin(phi), abs(Ef)); shading flat; colormap('jet');
% 
%   figure;
%   surface(rho .* cos(phi), rho .* sin(phi), abs(Hf));  shading flat; colormap('jet');

end
Corr_FEKO;

 figure(1);
 hold on;
 plot(rho(1, :), abs(E_tot_reshape(:, 1)));
 
 
Corr_FEKO_H;

 figure(2);
 hold on;
 plot(rho(1, :), abs(H_tot_reshape(:, 1)));
