c0 = 3e8;

F = 14e9; % Frequency at which far field is requested

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

rt = 0.0405319403216/1.9;
rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;
rd = 2.2e-2;
re = 2.3e-2;

R = [rr rp rt rd re]; % radius vector

drho = R(end)/100;
dphi = pi/180;

[rho, ph] = meshgrid(eps:drho:R(end), eps:dphi:2*pi+eps);

n = length(R);

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

epsilon = er .* er0;
mu = mur .* mu0;

L = 1e-3 * [1 1 1 20 0.5];
k = 10;

[E_aperture_rho, E_aperture_phi] = Near_field_fun_2(er, mur, R, F, L, rho, ph, k, drho, dphi);


%% Magnetic equivalent current on the aperture

[theta, phi] = meshgrid(-pi-eps:pi/180:pi+eps, eps:pi/180:2*pi+eps);

[Eth, Eph] = FF_M_V2(E_aperture_rho, E_aperture_phi, rho, ph, theta, phi, F, drho, dphi);

E_abs = sqrt(abs(Eth).^2 + abs(Eph).^2);

    

figure(41);

hold on;
plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))/max(abs(E_abs(1, :)))), 'LineWidth', 2);
hold on;
plot(theta(91, :)*(180/pi), db(abs(E_abs(91, :))/max(abs(E_abs(90, :)))), 'LineWidth', 2);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;

ylim([-50 0]);