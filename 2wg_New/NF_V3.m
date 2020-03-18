clear;

F = 14e9;

err = 1; erp = 1; murr = 1; murp = 1;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;       % Permeability in the medium
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;

rr = 0.0405319403216/2.1;
rp = 0.0405319403216/2;


[fc_1] = fc(rr, err, murr);
[fc_2] = fc(rp, erp, murp);

Nr1 = find(fc_1 < F);
Np1 = find(fc_2 < F);

Nr = 1:1:10;
Np = 1:1:10;

[X_] = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr); 

[Spp_, Spr_, Srp_, Srr_] = GSM(Nr, Np, F, rp, rr, erp, murp, err, murr, X_);

Slp = SL(rr, F, Np, 0.0005);
Slr = SL(rp, F, Nr, 0.0005);

Spp = Slp * Spp_ * Slp;
Spr = Slp * Spr_ * Slr;
Srp = Slr * Srp_ * Slp;
Srr = Slr * Srr_ * Slr;

%% 
[rho, phi] = meshgrid(eps:rr/100:rr,  eps:pi/180:2*pi-eps);

z = 0;

[Er_rho, Er_phi, Er_z] = E_r(Nr1, rho, phi, F, rr, z, epsilonr, mur);

ap = zeros(length(Np), 1);
ar = ones(length(Nr), 1);

br = Srp * ap + Srr * ar;

Gamma_sum = ar + br;

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));

for k = 1:length(Nr1)
    E_aperture_rho = E_aperture_rho + squeeze(Er_rho(k, :, :)) .* abs(Gamma_sum(k));
    E_aperture_phi = E_aperture_phi + squeeze(Er_phi(k, :, :)) .* abs(Gamma_sum(k));
    E_aperture_z = E_aperture_z + squeeze(Er_z(k, :, :)) .* abs(Gamma_sum(k));
end

% E_aperture = sqrt(squeeze(abs(Er_rho(5, :, :))).^2 + squeeze(abs(Er_phi(5, :, :))).^2 + squeeze(abs(Er_z(5, :, :))).^2);
E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2); % + abs(E_aperture_z).^2);

% % E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2); % + abs(E_aperture_z).^2);
x = rho .* cos(phi);
y = rho .* sin(phi);

figure;

surface(x,y, db(abs(E_aperture)./max(abs(E_aperture)))); shading flat;
colormap('jet');
figure;

surface(x,y, db(abs(E_aperture_rho)./max(abs(E_aperture_rho)))); shading flat;
colormap('jet');
figure;

surface(x,y, db(abs(E_aperture_phi)./max(abs(E_aperture_phi)))); shading flat;
colormap('jet');
figure;

surface(x,y, db(abs(E_aperture_z)./max(abs(E_aperture_z)))); shading flat;
colormap('jet');
plot_feko_NF;


%%  ----------------------------------------------------------------------------------------------------------

[rho, phi] = meshgrid(eps:rp/100:rp, eps:pi/180:2*pi+eps);

z = 0;

[Ep_rho, Ep_phi, Ep_z] = E_r(Np1, rho, phi, F, rp, z, epsilonp, mup);

ap = zeros(length(Np), 1);
ar = ones(length(Nr), 1);

bp = Spp * ap + Spr * ar;
Gamma_sum = ap + bp;

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));

for k = 1:length(Np1)
    E_aperture_rho = E_aperture_rho + squeeze(Ep_rho(k, :, :)) .* (Gamma_sum(k));
    E_aperture_phi = E_aperture_phi + squeeze(Ep_phi(k, :, :)) .* (Gamma_sum(k));
    E_aperture_z = E_aperture_z + squeeze(Ep_z(k, :, :)) .* (Gamma_sum(k));
end

% E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2 + abs(E_aperture_z).^2);

E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2);

x = rho .* cos(phi);
y = rho .* sin(phi);
% 
figure;
surface(x, y, db((abs(E_aperture))./max(abs(E_aperture)))); shading flat;
colormap('jet');
figure;

surface(x,y, db(abs(E_aperture_rho)./max(abs(E_aperture_rho)))); shading flat;
colormap('jet');
figure;
surface(x,y, db(abs(E_aperture_phi)./max(abs(E_aperture_phi)))); shading flat;
colormap('jet');
figure;

surface(x,y, db(abs(E_aperture_z)./max(abs(E_aperture_z)))); shading flat;
colormap('jet');
plot_feko_NF;
