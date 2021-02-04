function [yin, Gamma_K] = Mishustin_K_Freq_TE12(R, k0, L, er, mur, len, n_mode)

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

c0 = 3e8;

% beta = k0;
omega = c0.*k0;
% F = omega/(2 * pi);

% Str = load('Xmn_azimuthal_inc_TE.mat');
% 
% str = Str.xmn_TE;

% for i = 1:length(R)

r = R;
% 
Str = load('Xmn_azimuthal_inc_TE.mat');
Xmn = Str.xmn_TE;
xmn = Xmn(n_mode).xmn;

for k = 1:length(k0)

C1 = 2.*(xmn)^2.*(xmn./(k0(k) .* r))^2./((xmn.^2 - 1).*sqrt(1 - (xmn./(k0(k) .*r))^2));

C2 = 2./((xmn.^2 - 1).*sqrt(1 - (xmn./(k0(k) .* r)).^2));

I1 = @(beta) beta .* (-1j.*sqrt(-(1 - beta.^2))) .* (besselj_der(1, k0(k).*r.*beta)).^2./((xmn./(k0(k) .* r)).^2 - beta.^2).^2;
I2 = @(beta) (besselj(1, k0(k) .* r .* beta)).^2./(beta .* (-1j.*sqrt(-(1 - beta.^2))));

del = 0.01.*k0(k);

wpoints = [eps (1+1j).*eps k0(k) + 1j.*eps];

func = @(beta) C1 .* I1(beta) + C2 .* I2(beta);

% xmn11 = str(1).xmn;
% beta_rho11 = xmn11./r;
% 
% beta_z11 = -1j .* sqrt(-(k0(k).^2 - beta_rho11.^2));

% yin = conj(integral(func, eps, 50.*k0+1j.*eps, 'Waypoints', wpoints));

yin_(k) = integral(func, eps, 1+1j.*eps);

yin_1(k) = integral(func, 1+1j.*eps, L*k0(k)+1j.*eps);

yin(k) = real(yin_(k)) + 1j.*imag(yin_1(k));

% yap(i) = yin(i) .* sqrt(1 - (xmn./(k0.*r))^2);


% Gamma(i) = (1 - yap(i))/(1 + yap(i));
Gamma_K(k) = (1 - yin(k))./(1 + yin(k));

end

end