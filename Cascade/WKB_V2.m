
%% Internal impedance of TE11
clear;
c0 = 3e8;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

er = 1; mur = 1; epsilon = er .* er0; mu = mur .* mu0;

str = load('Xmn.mat');
str = str.Xmn;

% str = load('Xmn_azimuthal_inc_TE.mat');
% str = str.xmn_TE;

F = 5e9;
lamb = c0./F;

beta = 2.*pi./lamb;

omega = 2 .* pi .* F;


rr = 2e-2;
rt = 2.*lamb;

d = 6 .* lamb;

% slope = (r(1) - r(end))./d;


fc_ = fc(rt, er, mur);

fc_end = find(fc_ < F);

N = length(fc_end);

for i = 1:N
    fc_i = fc_(i);
    r(i) = str(N - i + 1).xmn .*c0./(2 .* pi .* F);
    Num(i) = i;
end

r_vec = [rt r];

slope = (r_vec(1) - r_vec(end))./d;

z_axis = linspace(eps, d, length(r_vec));

Psi_0 = acos(str(N).xmn./(beta .* r_vec(1)));


drho = r(1)/100;
dphi = pi./180;
phi = linspace(eps, 2 .* pi, 360);

[rho, phi] = meshgrid(eps:drho:r_vec(1), phi);
x = rho .* cos(phi);
y = rho .* sin(phi);

dth = pi/180;
dph = pi/180;

[theta_obs, phi_obs] = meshgrid(-pi/2-eps:dth:pi/2-eps, eps:dph:2*pi+eps);

Erho = zeros(size(rho));
Ephi = zeros(size(rho));


Eth = zeros(size(theta_obs));
Eph = zeros(size(theta_obs));
Eco = zeros(size(theta_obs));
Exp = zeros(size(theta_obs));
% u = find(r(Num == 10));

for k = 1:N

    
m = str(N - k + 1).m;

Psi_z = acos(str(N - k + 1).xmn./(beta .* r_vec(N - k + 1)));

beta_rho = str(N - k + 1).xmn./r_vec(1);
beta_z = -1j .* sqrt(-(beta.^2 - (beta_rho).^2));

h = (2./sqrt(beta_z)) .* exp(-1j .* str(N - k + 1).xmn .*((tan(Psi_z) - Psi_z)) - (tan(Psi_0) - Psi_0));
% h = 1;

f = besselj(m, beta_rho * rho);
g = cos(m .* phi);

Pot_00 = h .* f .* g;


if str(N - k + 1).mode == "TM"
    
  K_00 = 1j .* beta_z./(beta_rho).^2;

  Erho_i = -K_00 .* beta_rho .* besselj_der(m, beta_rho .* rho) .* cos(m .* phi) .* h;
   
  Ephi_i = K_00 .* m./rho .* besselj(m, beta_rho .* rho) .* sin(m .* phi) .* h;
  
  [Eth_i, Eph_i, Eco_i, Exp_i] = FF_M_V2(Erho_i, Ephi_i, rho, phi, theta_obs, phi_obs, F, drho, dphi);
 
  
else
   Z_00 = (omega .* mu)./beta_z;
   
   K_00 = 1j .* beta_z./(beta_rho).^2 .* Z_00;

   Ephi_i = K_00 .* beta_rho .* besselj_der(m, beta_rho .* rho) .* cos(m .* phi) .* h;
   
   Erho_i = K_00 .* m./rho .* besselj(m, beta_rho .* rho) .* sin(m .* phi) .* h;
   
   if m == 1
   
      [Eth_i, Eph_i, Eco_i, Exp_i, CO_i, XP_i] = FF_apertureFSCir3(str(N - k + 1).n, 1, 0, theta_obs, phi_obs, F, er, mur, r_vec(1));
   
   else
      [Eth_i, Eph_i, Eco_i, Exp_i] = FF_M_V2(Erho_i, Ephi_i, rho, phi, theta_obs, phi_obs, F, drho, dphi);
   end
    
end

Erho = Erho + Erho_i;
Ephi = Ephi + Ephi_i;

Eth = Eth + h .* Eth_i;
Eph = Eph + h .* Eph_i;

Eco = Eco + h .* Eco_i;
Exp = Exp + h .* Exp_i;


end



E_aper = sqrt(abs(Erho).^2 + abs(Ephi).^2);
E_FF = sqrt(abs(Eth).^2 + abs(Eph).^2);

% figure;
% 
% surface(x, y, abs(Pot_00)); shading flat; colormap('jet');
figure;

surface(x, y, db(abs(Erho)./max(max(abs(Erho))))); shading flat; colormap('jet');

figure;

surface(x, y, db(abs(Ephi)./max(max(abs(Ephi))))); shading flat; colormap('jet');


figure;

surface(x, y, db(abs(E_aper)./max(max(abs(E_aper))))); shading flat; colormap('jet');

figure;
plot(theta_obs(1, :)*180./pi, db(E_FF(1, :)./max(((abs(E_FF(1, :))))))); hold on;

plot(theta_obs(91, :).*180/pi, db(E_FF(91, :)./max(((abs(E_FF(91, :))))))); hold on; grid on; ylim([-50 0]);

figure;
plot(theta_obs(1, :)*180./pi, db(Eco(1, :)./max(((abs(Eco(1, :))))))); hold on;

plot(theta_obs(46, :).*180/pi, db(Exp(46, :)./max(((abs(Eco(1, :))))))); hold on; grid on; ylim([-50 0]);



