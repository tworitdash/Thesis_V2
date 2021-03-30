
f_req = 5e9;
er = 1;
mur = 1;
N = 1;
r_ = 1.8e-2;
[rho, phi] = meshgrid(linspace(eps, r_, 100), linspace(eps, 2*pi, 360));

drho = r_/100;
dphi = pi/180;

[erho_, ephi_] = E_n(1:1:N(end), rho, phi, f_req, R(end), z, epsilon(end), mu(end), drho, dphi);
erho = squeeze(erho_(1, :, :));
ephi = squeeze(ephi_(1, :, :));
eres = sqrt(abs(erho).^2 + abs(ephi).^2);

x = rho .* cos(phi);
y = rho .* sin(phi);
%% 
figure;

surface(x, y, (abs(erho))); shading flat;

colormap('jet');

figure;

surface(x, y, (abs(ephi))); shading flat;

colormap('jet');

figure;

surface(x, y, (abs(eres))); shading flat;

colormap('jet');