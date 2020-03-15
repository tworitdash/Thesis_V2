clear;

rr = 2e-2; % Base radius
rt = 4e-2; % Top redius
n = 20; % number of transitions


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

R = linspace(rr, rt, n); % radius vector

n = length(R);

F = 14e9; % Frequency of operation

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG section

epsilon = er .* er0;
mu = mur .* mu0;


Length = 5e-2; % height of the cone

L = ones(1, n) .* Length/n; % length of each waveguide section
L(1) = L(1)/2;
L(end) = L(end)/2;

[STT, STR, SRT, SRR, N] = GSM_N(R, L, er, mur, F);

%% 
[rho, phi] = meshgrid(eps:R(1)/100:R(1),  eps:pi/180:2*pi-eps);

z = 0;

[Er_rho, Er_phi, Er_z] = E_r(1:1:N(1), rho, phi, F, R(1), z, epsilon(1), mu(1));

ap = zeros(N(end), 1);
ar = ones(N(1), 1);


br = squeeze(SRT(1, :, :)) * ap + squeeze(SRR(1, :, :)) * ar;

Gamma_sum = ar + br;

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));

for k = 1:N(1)
    E_aperture_rho = E_aperture_rho + squeeze(Er_rho(k, :, :)) .* abs(Gamma_sum(k));
    E_aperture_phi = E_aperture_phi + squeeze(Er_phi(k, :, :)) .* abs(Gamma_sum(k));
    E_aperture_z = E_aperture_z + squeeze(Er_z(k, :, :));% .* abs(Gamma_sum(k));
end

% E_aperture = sqrt(squeeze(abs(Er_rho(5, :, :))).^2 + squeeze(abs(Er_phi(5, :, :))).^2 + squeeze(abs(Er_z(5, :, :))).^2);
E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2 + abs(E_aperture_z).^2);
% E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2); % + abs(E_aperture_z).^2);


x = rho .* cos(phi);
y = rho .* sin(phi);

% figure;
% 
% surface(x,y, db(abs(E_aperture)./max(E_aperture))); shading flat;
% 
% figure;
% 
% surface(x,y, db(abs(E_aperture_rho)./max(E_aperture_rho))); shading flat;
% 
% figure;
% 
% surface(x,y, db(abs(E_aperture_phi)./max(E_aperture_phi))); shading flat;
% 
figure;

surface(x,y, db(abs(E_aperture_z)./max(E_aperture_z))); shading flat;

colormap('jet');

plot_feko_NF;


%%  ----------------------------------------------------------------------------------------------------------

[rho, phi] = meshgrid(eps:R(end)/100:R(end), eps:pi/180:2*pi+eps);

z = 0;

[Ep_rho, Ep_phi, Ep_z] = E_r(1:1:N(end), rho, phi, F, R(end), z, epsilon(end), mu(end));


ap = zeros(N(end), 1);
ar = ones(N(1), 1);

bp = squeeze(STT(1, :, :)) * ap + squeeze(STR(1, :, :)) * ar;
Gamma_sum = ap + bp;

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));

for k = 1:N(end)
    E_aperture_rho = E_aperture_rho + squeeze(Ep_rho(k, :, :)) .* (Gamma_sum(k));
    E_aperture_phi = E_aperture_phi + squeeze(Ep_phi(k, :, :)) .* (Gamma_sum(k));
    E_aperture_z = E_aperture_z + squeeze(Ep_z(k, :, :)) .* (Gamma_sum(k));
end

E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2 + abs(E_aperture_z).^2);

% E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2);

x = rho .* cos(phi);
y = rho .* sin(phi);
% 
figure;
surface(x, y, db((abs(E_aperture))./max(abs(E_aperture)))); shading flat;

figure;

surface(x,y, db(abs(E_aperture_rho)./max(E_aperture_rho)));


figure;

surface(x,y, db(abs(E_aperture_phi)./max(E_aperture_phi)));


figure;

surface(x,y, db(abs(E_aperture_z)./max(E_aperture_z)));

plot_feko_NF;





