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

L = 1e-3 * [1 1 1 20 1];

[E_aperture_rho, E_aperture_phi, E_aperture_z] = Near_field_fun(er, mur, R, F, L, rho, ph);


%% Magnetic equivalent current on the aperture

M_rho = E_aperture_phi;
M_phi = -E_aperture_rho;

Mx = cos(ph) .* M_rho - sin(ph) .* M_phi;
My = -sin(ph) .* M_rho + cos(ph) .* M_phi;
Mz = zeros(size(rho));

[theta, phi] = meshgrid(-pi+eps:pi/180:pi+eps, eps:pi/180:2*pi+eps);
    
k0 = (2 * pi * F)/c0;

kx = k0 * sin(theta) .* cos(phi);
ky = k0 * sin(theta) .* sin(phi);

kz = (-1j)*sqrt(-(k0.^2 - kx.^2 - ky.^2));
    
zeta = 120*pi;

c_new = (-1./(2 * zeta * k0 * kz));

[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = ...
    SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c_new);

[M_ft_x, M_ft_y, M_ft_z] = J_Aperture(Mx, My, Mz, drho, dphi, k0, rho, ph, theta, phi);

[Hx, Hy, Hz] = FF(c_new, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, M_ft_x, M_ft_y, M_ft_z);

H_abs = sqrt(Hx.^2 + Hy.^2 + Hz.^2);

figure;

    
plot(theta(1, :)*(180/pi), db(abs(H_abs(1, :))/max(abs(H_abs(1, :)))), 'LineWidth', 2);
 
xlabel('\theta(Deg), \theta_0(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs}(\phi = 0) (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far magnetic field', 'FontSize', 12, 'FontWeight', 'bold');
            
grid on;

