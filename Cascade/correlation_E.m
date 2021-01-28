c0 = 3e8;


er = 1;
mur = 1;
R = 1.8e-2;
F = 17e9;
lamb = c0/F; 

drho = R/100;
dphi = 180/pi;

L = lamb/4;
[rho, phi] = meshgrid(eps:R/100:R, eps:pi/180:2*pi+eps);

x = rho .* cos(phi);
y = rho .* sin(phi);

N_modes = 9:1:10;


str = load('Xmn.mat');

str = str.Xmn;



% [E_rho, E_phi] = Near_field_fun_2(er, mur, R, F, L, rho, phi, 0, drho, dphi);

[E_rhoTheory, E_phiTheory] = E_n(N_modes, rho, phi, F, R(end), L, epsilon(end), mu(end), drho, dphi);
%% Near fidlds from FEKO open waveguide

% excitation
Corr_FEKO;

figure; surface(x_f, y_f, ((abs(E_tot_reshape).')./max(max(abs(E_tot_reshape).')))); shading flat; colormap('jet');


for i = 1:length(N_modes)
    mode = str(N_modes(i)).mode;
    m = str(N_modes(i)).m;
    n = str(N_modes(i)).n;
    E_rho_i = squeeze(E_rhoTheory(i, :, :));
    E_phi_i = squeeze(E_phiTheory(i, :, :));
    E_i = abs(sqrt(abs(E_rho_i).^2 + abs(E_phi_i).^2));
    figure; surface(x, y, (abs(E_i)./(max(max(abs(E_i)))))); shading flat; colorbar; colormap('jet'); title(['E_{abs} of ', mode, '_{', num2str(m), num2str(n), '}'])
    
%% Correlation

    Corr = E_i(1:end-1, :) .* E_tot_reshape.';

    figure; surface(x_f, y_f, abs(Corr)./max(max(Corr)));shading flat;colormap('jet')
%     figure; surface(rho(1, :), phi(:, 1)*180/pi, db(abs(E_rho))); shading flat; colorbar('jet'); 
%     figure; surface(rho(1, :), phi(:, 1)*180/pi, db(abs(E_phi))); shading flat; colorbar('jet');
end

 

