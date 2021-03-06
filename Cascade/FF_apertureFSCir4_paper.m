function [Eth, Eph, Eco, Exp, CO, XP] = FF_apertureFSCir4_paper(N, Dm, Gamma, theta, phi, F, er, mur, R, R1) 

c0 = 3e8;
v = c0./sqrt(er);

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

beta = 2 * pi * F / v;


Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;

Eth = zeros(size(theta));
Eph = zeros(size(theta));

Eco = zeros(size(theta));
Exp = zeros(size(theta));

CO = zeros(size(theta));
XP = zeros(size(theta));

for i = 1:length(N)

xmn = str(N(i)).xmn;
m = str(N(i)).m;
mode = str(N(i)).mode;

% beta_rhop = xmn./R;

beta_rhop = xmn./R;
% beta_rhop1 = xmn./R1;
beta_z = -1j.*sqrt(-(beta.^2 - beta_rhop.^2));

ZTE = 2 * pi * F * mu ./ beta_z;

Nup = 1j .* beta_z ./ beta_rhop.^2 * ZTE;
% Nup = 1;
k0 = (2 * pi * F)/c0;

kx = k0 .* sin(theta) .* cos(phi);
ky = k0 .* sin(theta) .* sin(phi);

kz = (-1j)*sqrt(-(k0.^2 - kx.^2 - ky.^2));

Nu = sqrt((kx.^2 + ky.^2));

% Psi = acos(kx./Nu);
Psi = phi;

for l = 1:size(Nu, 1)
    
    for k = 1:size(Nu, 2)

        I0(l, k) = Lommel(0, R, beta_rhop, Nu(l, k), m - 1, m - 1);
        I2(l, k) = Lommel(0, R, beta_rhop, Nu(l, k), m + 1, m + 1);

    end
end

E_ft_x = (Nup) .* Dm(i) .* pi .* beta_rhop .* sin(2 * Psi) .* I2;
E_ft_y = (Nup) .* Dm(i) .* pi .* beta_rhop .* (I0 + cos(2 .* Psi) .* I2);


c2 = 1j .* k0 / (4 * pi);

Eth_o = c2 .* (1 + cos(theta)) .* (E_ft_x .* cos(phi) + E_ft_y .* sin(phi)) .* (1 + Gamma);
Eph_o = c2 .* (1 + cos(theta)) .* (E_ft_y .* cos(phi) - E_ft_x .* sin(phi)) .* (1 + Gamma);


Eth = Eth + c2 .* (1 + cos(theta)) .* (E_ft_x .* cos(phi) + E_ft_y .* sin(phi)) .* (1 + Gamma);
Eph = Eph + c2 .* (1 + cos(theta)) .* (E_ft_y .* cos(phi) - E_ft_x .* sin(phi)) .* (1 + Gamma);

Eres = sqrt(abs(Eth_o).^2 + abs(Eph_o).^2);
% figure; plot(theta(1, :)*180/pi, db(abs(Eres)));
% Exp = Exp +  (1 + cos(theta)) .* E_ft_x;
% Eco = Eco +  (1 + cos(theta)) .* E_ft_y;

Exp = Exp + cos(phi) .* Eth_o - sin(phi) .* Eph_o;
Eco = Eco + sin(phi) .* Eth_o + cos(phi) .* Eph_o;

XP = XP + (Nup) .* Dm(i) .* pi .* beta_rhop .* I2 .* cos(theta);
CO = CO + (Nup) .* Dm(i) .* pi .* beta_rhop .* I0 .* cos(theta);


end

end