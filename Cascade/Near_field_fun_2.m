function [E_aperture_rho, E_aperture_phi] = Near_field_fun_2(er, mur, R, F, L, rho, phi, k, drho, dphi)


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er .* er0;
mu = mur .* mu0;

[STT, STR, SRT, SRR, N] = GSM_N(R, L, er, mur, F, k);


%%  ----------------------------------------------------------------------------------------------------------

z = 0;

[Ep_rho, Ep_phi] = E_n(1:1:N(end), rho, phi, F, R(end), z, epsilon(end), mu(end), drho, dphi);

% [Ep_rho, Ep_phi, ~] = E_r(1:1:N(end), rho, phi, F, R(end), z, epsilon(end), mu(end));


ap = zeros(N(end), 1);
ar = ones(N(1), 1);

if N(1) == 1
    STR_req = STR(1, 1:N(end));
else
    STR_req = squeeze(STR(1, 1:N(end), 1:N(1)));
end

bp = squeeze(STT(1, 1:N(end), 1:N(end))) * ap + STR_req * ar;
Gamma_sum = ap + bp;

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
% E_aperture_z = zeros(size(rho));

for k = 1:N(end)
    
    E_aperture_rho = E_aperture_rho + squeeze(Ep_rho(k, :, :)) .* (Gamma_sum(k)); %.* exp(1j .* pi/2);
    E_aperture_phi = E_aperture_phi + squeeze(Ep_phi(k, :, :)) .* (Gamma_sum(k)); %.* exp(1j .* pi/2);
%     E_aperture_z = E_aperture_z + squeeze(Ep_z(k, :, :)) .* (Gamma_sum(k));
end

% E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2 + abs(E_aperture_z).^2);

% E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2);

% x = rho .* cos(phi);
% y = rho .* sin(phi);
% 
% figure;
% surface(x, y, db((abs(E_aperture))./max(max(abs(E_aperture))))); shading flat;
% 
% % surface(x, y, (abs(E_aperture))); shading flat;
% 
% colormap('jet');
% figure;
% 
% surface(x,y, db(abs(E_aperture_rho)./max(max(abs(E_aperture_rho))))); shading flat;
% colormap('jet');
% 
% figure;
% 
% surface(x,y, db(abs(E_aperture_phi)./(max(max(abs(E_aperture_phi)))))); shading flat;
% colormap('jet');

end