
clear;
c0 = 3e8;


er = 1;
mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

R = 1.8e-2;
F = 17e9;
lamb = c0/F; 

drho = R/100;
dphi = 180/pi;

L = lamb/4;
[rho, phi] = meshgrid(1e-7:R/100:R, linspace(1e-7, 2*pi, 360));

x = rho .* cos(phi);
y = rho .* sin(phi);

N_modes = [9 10];


str = load('Xmn.mat');

str = str.Xmn;


[E_rhoTheory, E_phiTheory] = E_n(N_modes, rho, phi, F, R(end), L, epsilon(end), mu(end), drho, dphi);
%% Near fidlds from FEKO open waveguide

% excitation
Corr_FEKO;

% figure; surface(x_f, y_f, (abs(E_tot_reshape).')./max(max(abs(E_tot_reshape).'))); shading flat; colormap('jet');

figure; surface(x_f, y_f, (abs(E_tot_reshape).')); shading flat; colormap('jet'); colorbar;

% figure; surface(x_f, y_f, (abs(E_rho_reshape).')); shading flat; colormap('jet'); colorbar;


% figure; surface(x_f, y_f, (abs(E_phi_reshape).')); shading flat; colormap('jet'); colorbar;

 for i = 1:length(N_modes)
     mode = str(N_modes(i)).mode;
     m = str(N_modes(i)).m;
     n = str(N_modes(i)).n;
     E_rho_i = squeeze(E_rhoTheory(i, :, :));
     E_phi_i = squeeze(E_phiTheory(i, :, :));
     E_i = abs(sqrt(abs(E_rho_i).^2 + abs(E_phi_i).^2));
%      figure; surface(x, y, (abs(E_i)./(max(max(abs(E_i)))))); shading flat; colorbar; colormap('jet'); title(['E_{abs} of ', mode, '_{', num2str(m), num2str(n), '}'])
     figure; surface(x, y, (abs(E_i))); shading flat; colorbar; colormap('jet'); title(['E_{abs} of ', mode, '_{', num2str(m), num2str(n), '}'])
%      figure; surface(x, y, (abs(E_rho_i))); shading flat; colorbar; colormap('jet'); title(['E_{abs} of ', mode, '_{', num2str(m), num2str(n), '}'])
%      figure; surface(x, y, (abs(E_phi_i))); shading flat; colorbar; colormap('jet'); title(['E_{abs} of ', mode, '_{', num2str(m), num2str(n), '}'])
end


b = (squeeze(E_phiTheory(1, :, :)) .* E_rho_reshape.' - squeeze(E_rhoTheory(1, :, :)) .* E_phi_reshape.')./...
    (squeeze(E_phiTheory(1, :, :)) .* squeeze(E_rhoTheory(2, :, :)) - squeeze(E_rhoTheory(1, :, :)) .* squeeze(E_phiTheory(2, :, :)));
 
% a = (E_rho_reshape.' - b .* squeeze(E_rhoTheory(2, :, :)))./squeeze(E_rhoTheory(1, :, :));
 
a = (E_rho_reshape.')./squeeze(E_rhoTheory(1, :, :));

% figure; surface(x, y, abs(b./abs(max(max(b))))); shading flat;colorbar; colormap('jet');

% figure; surface(x, y, abs(b)); shading flat;colorbar; colormap('jet');
% 
% % figure; surface(x, y, abs(a./abs(max(max(a))))); shading flat;colorbar; colormap('jet');
% 
% figure; surface(x, y, abs(a)); shading flat;colorbar; colormap('jet');