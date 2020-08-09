function [Eth, Eph] = FF_M_cart_numerical(Eax, Eay, Eaz, x, y, z, theta, phi, F, dx, dy)

c0 = 3e8;


    
k0 = (2 * pi * F)/c0;

% kx = k0 .* sin(theta) .* cos(phi);
% ky = k0 .* sin(theta) .* sin(phi);
% 
% kz = (-1j)*sqrt(-(k0.^2 - kx.^2 - ky.^2));
    
% zeta = 120*pi;
% 
% c_new = (-1./(2 * zeta * k0 * kz));
% 
% [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);
% 
% [SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = ...
%     SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c_new);

[E_ft_x, E_ft_y, E_ft_z] = J_cart(Eax, Eay, Eaz, dx, dy, k0, x, y, z, theta, phi);
% [E_ft_x, E_ft_y, E_ft_z] = J_cart_analytical(E0, k0, theta, phi, a, b, M, Dm, Gamma);

r_obs = 10000 * 2 * pi * F ./ c0;
% 
% c2 = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

% [Hx, Hy, Hz] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, M_ft_x, M_ft_y, M_ft_z);

% c2 = 1j .* k0 .* exp(-1j * k0 * r_obs) / (4 * pi);
% 
% Eth = c2 .* sin(phi) .* (E_ft_x + E_ft_y + E_ft_z);
% Eph = c2 .* cos(theta) .* cos(phi) .* (E_ft_x + E_ft_y + E_ft_z);

c2 = 1j .* k0 .* exp(-1j * k0 * r_obs) / (4 * pi);
% 
% Eth = c2 .* (E_ft_x .* cos(phi) + E_ft_y .* sin(phi));
% Eph = c2 .* cos(theta) .* (E_ft_y .* cos(phi) - E_ft_x .* sin(phi));


Eth = c2 .* sin(phi) .* (E_ft_y);
Eph = c2 .* cos(theta) .* cos(phi) .* (E_ft_y);

% Ex = zeta * Hy;
% Ey = -zeta * Hx;
% 
% Ez = 0;


end