function [Ex, Ey, Ez, Hx, Hy, Hz] = FF_M(E_aperture_rho, E_aperture_phi, rho, ph, theta, phi, F, drho, dphi)

c0 = 3e8;

M_rho = E_aperture_phi;
M_phi = -E_aperture_rho;

Mx = cos(ph) .* M_rho - sin(ph) .* M_phi;
My = sin(ph) .* M_rho + cos(ph) .* M_phi;
Mz = zeros(size(rho));


    
k0 = (2 * pi * F)/c0;

kx = k0 .* sin(theta) .* cos(phi);
ky = k0 .* sin(theta) .* sin(phi);

kz = (-1j)*sqrt(-(k0.^2 - kx.^2 - ky.^2));
    
zeta = 120*pi;

c_new = (-1./(2 * zeta * k0 * kz));

[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = ...
    SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c_new);

[M_ft_x, M_ft_y, M_ft_z] = J_Aperture(Mx, My, Mz, drho, dphi, k0, rho, ph, theta, phi);

r_obs = 10000 * 2 * pi * F ./ c0;

c2 = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

[Hx, Hy, Hz] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, M_ft_x, M_ft_y, M_ft_z);


Ex = zeta * Hy;
Ey = -zeta * Hx;

Ez = 0;


end