
%% 
% clear;
c0 = 3e8;
er = 1;
mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;



N_modes = 9:100;
R = 1.8e-2;
F = 17e9;

lamb = c0/F;
L = lamb/5;

[rho, phi] = meshgrid(linspace(eps, R, 100), linspace(eps, 2 * pi, 360));

drho = rho(1, 2) - rho(1, 1);
dphi = phi(2, 1) - phi(1, 1);

[E_rho, E_phi, ~, H_rho, H_phi] = E_n(N_modes, rho, phi, F, R(end), L, epsilon(end), mu(end), drho, dphi);
Corr_FEKO;
Corr_FEKO_H;
Power_F = zeros(1, length(N_modes));

Power_F_rev = zeros(1, length(N_modes));


for j = 1:length(N_modes)
   H_rho_j = squeeze(H_rho(j, :, :));
   H_phi_j = squeeze(H_phi(j, :, :));
   E_rho_j = squeeze(E_rho(j, :, :));
   E_phi_j = squeeze(E_phi(j, :, :));
   Power_F(j) = Pow_2modes(E_rho_reshape_.', E_phi_reshape_.', H_rho_j, H_phi_j, rho_f, drho, dphi);
   Power_F_rev(j) = Pow_2modes(E_rho_j, E_phi_j, H_rho_reshape_.', H_phi_reshape_.', rho_f, drho, dphi);
end

diff = (Power_F - Power(9, 9:end));

figure; plot(9:100, db(abs(diff))/2);

figure; plot(9:100, db(abs(Power(9, 9:100)))/2); hold on; plot(9:100, db(abs(Power_F)));


figure; plot(9:100, db(abs(Power_F./max(Power_F))));

% Power_F = Pow_2modes(E_rho_reshape_.', E_phi_reshape_.', H_rho_reshape_.', H_phi_reshape_.', rho_f, drho, dphi);
