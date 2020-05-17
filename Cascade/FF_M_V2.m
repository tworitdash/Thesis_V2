function [Eth, Eph] = FF_M_V2(E_aperture_rho, E_aperture_phi, rho, ph, theta, phi, F, drho, dphi)

c0 = 3e8;



Eax = cos(ph) .* E_aperture_rho - sin(ph) .* E_aperture_phi;
Eay = sin(ph) .* E_aperture_rho + cos(ph) .* E_aperture_phi;
Eaz = zeros(size(rho));


    
k0 = (2 * pi * F)/c0;

kx = k0 .* sin(theta) .* cos(phi);
ky = k0 .* sin(theta) .* sin(phi);

kz = (-1j)*sqrt(-(k0.^2 - kx.^2 - ky.^2));
    
% zeta = 120*pi;
% 
% c_new = (-1./(2 * zeta * k0 * kz));
% 
% [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);
% 
% [SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = ...
%     SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c_new);

[E_ft_x, E_ft_y, E_ft_z] = J_Aperture(Eax, Eay, Eaz, drho, dphi, k0, rho, ph, theta, phi);

r_obs = 10000 * 2 * pi * F ./ c0;
% 
% c2 = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

% [Hx, Hy, Hz] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, M_ft_x, M_ft_y, M_ft_z);

c2 = 1j .* k0 .* exp(-1j * k0 * r_obs) / (4 * pi);

Eth = c2 .* (E_ft_x .* cos(phi) + E_ft_y .* sin(phi));
Eph = c2 .* cos(theta) .* (E_ft_y .* cos(phi) - E_ft_x .* sin(phi));

% Ex = zeta * Hy;
% Ey = -zeta * Hx;
% 
% Ez = 0;


end