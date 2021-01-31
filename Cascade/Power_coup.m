
%% 
c0 = 3e8;
er = 1;
mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;



N_modes = 1:100;
R = 1.8e-2;
F = 17e9;

lamb = c0/F;
L = lamb/5;

[rho, phi] = meshgrid(linspace(eps, R, 100), linspace(eps, 2 * pi, 360));

drho = rho(1, 2) - rho(1, 1);
dphi = phi(2, 1) - phi(1, 1);

[E_rho, E_phi, ~, H_rho, H_phi] = E_n(N_modes, rho, phi, F, R(end), L, epsilon(end), mu(end), drho, dphi);
Power = zeros(length(N_modes), length(N_modes));
for i = 1:length(N_modes)
    E_rho_i = squeeze(E_rho(i, :, :));
    E_phi_i = squeeze(E_phi(i, :, :));
    for j = 1:length(N_modes)
        H_rho_j = squeeze(H_rho(j, :, :));
        H_phi_j = squeeze(H_phi(j, :, :));
        Power(i, j) = Pow_2modes(E_rho_i, E_phi_i, H_rho_j, H_phi_j, rho, drho, dphi);
    end
end

figure; surface(N_modes, N_modes, abs(Power));
