%% Far fields of the ctlindrical waveguide analytically for TE11 mode
clear;
c0 = 3e8;

er = 1; mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;


% 
% dth = pi/180;
% dphi = pi/180;

F = 5e9;

beta = 2 * pi * F ./ c0;

omega = 2 * pi * F;


R = 2e-2;

Str = load('Xmn.mat');

str = Str.Xmn;

xmn = str(9).xmn;
m = str(9).m;
mode = str(9).mode;

xmn1 = str(1).xmn;
m1 = str(1).m;
mode1 = str(1).mode;

beta_rhop = xmn./R;
beta_z = -1j .* sqrt(-(beta.^2 - beta_rhop.^2));

beta_rhop1 = xmn1./R;
beta_z1 = -1j .* sqrt(-(beta.^2 - beta_rhop1.^2));

ZTE = 2 * pi * F * mu ./ beta_z;
ZTE1 = 2 * pi * F * mu ./ beta_z1;

if m == 0
  deltam = 1;
else
  deltam = 0;
end
if mode == "TE"
   Nup = (pi*(1+deltam)/2 .* (xmn.^2 - m.^2) .* (besselj(m, xmn)).^2).^(-1);
   Nup1 = (pi*(1+deltam)/2 .* (xmn1.^2 - m1.^2) .* (besselj(m1, xmn1)).^2).^(-1);
elseif modep  == "TM"
   Nup = (pi*(1+deltam)/2 .* (xmn_p).^2 .* (besselj_der(m, xmn)).^2).^(-1);
end
% Nup = 1j .* beta ./ beta_rhop.^2 * ZTE;

k0 = (2 * pi * F)/c0;

% kx = k0 .* sin(theta) .* cos(phi);
% ky = k0 .* sin(theta) .* sin(phi);

kx = linspace(-5*k0, 5*k0, 1000);
ky = linspace(-5*k0, 5*k0, 1000);



[kx, ky] = meshgrid(kx, ky);

dkx = kx(1, 1) - kx(2, 1);
dky = ky(1, 1) - ky(1, 2);

Nu = sqrt(kx.^2 + ky.^2);


kz = (-1j)*sqrt(-(k0.^2 - Nu.^2));

% Nu = sqrt((kx.^2 + ky.^2));

Psi = acos(kx./Nu);

for i = 1:size(Nu, 1)
    
   for k = 1:size(Nu, 2)

        I0(i, k) = Lommel(0, R, beta_rhop, Nu(i, k), m - 1, m - 1);
        I2(i, k) = Lommel(0, R, beta_rhop, Nu(i, k), m + 1, m + 1);
        
        I0_2(i, k) = Lommel(0, R, beta_rhop1, Nu(i, k), m - 1, m - 1);
        I2_2(i, k) = Lommel(0, R, beta_rhop1, Nu(i, k), m + 1, m + 1);
   end
   
end

% E_ft_x = (Nup) .* pi .* beta_rhop .* sin(2 * Psi) .* I2;
% E_ft_y = (Nup) .* pi .* beta_rhop .* (I0 + cos(2 .* Psi) .* I2);
% 

E_ft_x =  pi .* beta_rhop .* sin(2 * Psi) .* I2;
E_ft_y =  pi .* beta_rhop .* (I0 + cos(2 .* Psi) .* I2);

E_ft_x1 = pi .* beta_rhop1 .* sin(2 * Psi) .* I2_2;
E_ft_y1 = pi .* beta_rhop1 .* (I0_2 + cos(2 .* Psi) .* I2_2);

Ixx = (k0^2 - ky^2)./(omega .* mu .* kz .* 4 .* pi^2) .* E_ft_x .* E_ft_x1;
Iyy = (k0^2 - kx^2)./(omega .* mu .* kz .* 4 .* pi^2) .* E_ft_y .* E_ft_y1;

Ixy = 2 .* kx .* ky ./ (omega .* mu .* kz .* 4 .* pi^2) .* E_ft_x .* E_ft_y1;

Exx = sum(sum(Ixx)) .* dkx .* dky;
Eyy = sum(sum(Iyy)) .* dkx .* dky;
Exy = sum(sum(Ixy)) .* dkx .* dky;

figure;
surf(kx(1, :)./k0, ky(:, 1)./k0, abs(Ixx), 'LineWidth', 2); shading flat;

figure;
surf(kx(1, :)./k0, ky(:, 1)./k0, abs(Iyy), 'LineWidth', 2); shading flat;

figure;
surf(kx(1, :)./k0, ky(:, 1)./k0, abs(Ixy), 'LineWidth', 2); shading flat;

