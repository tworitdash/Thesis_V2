
%% 
% clear;
c0 = 3e8;
er = 1;
mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;



N_modes = 1:100;
R = 3.*1.8e-2;
F = 16e9;

lamb = c0/F;
% L = lamb/5;
L = c0/5e9/4;

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

% diff = (Power_F - Power(1, :));

% figure; plot(:100, db(abs(diff))/2);

figure; plot(1:100, db(abs(Power(1, 1:100)))); hold on; plot(1:100, db(abs(Power_F)));


figure; plot(1:100, db(abs(Power_F./max(Power))));


grid on;
xlabel('Mode Number', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Cross Eergy in dB', 'FontWeight', 'bold', 'FontSize', 16);
title('Cross Energy. Integral of conjugate of Poynting vector', 'FontWeight', 'bold', 'FontSize', 16);

% Power_F = Pow_2modes(E_rho_reshape_.', E_phi_reshape_.', H_rho_reshape_.', H_phi_reshape_.', rho_f, drho, dphi);
