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

xmn = str(1).xmn;
m = str(1).m;
mode = str(1).mode;

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
elseif mode  == "TM"
   Nup = (pi*(1+deltam)/2 .* (xmn_p).^2 .* (besselj_der(m, xmn)).^2).^(-1);
end

k0 = (2 * pi * F)/c0;
%%

% Nu = linspace(0, 50.*beta_rhop1, 51);

% Nu = beta_rhop1;
% dnu = Nu(2) - Nu(1);

kz = @(Nu) -1j .* sqrt(-(k0.^2 - Nu.^2));

% for i = 1:length(Nu)

        I0 = @(Nu) Lommel2_h(0, R, beta_rhop, Nu, m - 1, m - 1);
        I2 =  @(Nu) Lommel2_h(0, R, beta_rhop, Nu, m + 1, m + 1);
        
        I0_2 = @(Nu) Lommel2_h(0, R, beta_rhop1, Nu, m - 1, m - 1);
        I2_2 = @(Nu) Lommel2_h(0, R, beta_rhop1, Nu, m + 1, m + 1);
   
% end

% Ixx = @(Nu) beta_rhop .* beta_rhop1 .* (1/(4) .* omega .* mu) .* pi/2 .* (2.*k0^2 - Nu.^2)/kz(Nu) .* I2(Nu) .* I2_2(Nu) .* Nu;
% Iyy = @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* (1/(4) .* omega .* mu) .* pi/2 .* (2.*k0^2 - Nu.^2)/kz(Nu) .* I0(Nu) .* I0_2(Nu) .* Nu;

I =  @(Nu) beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((2 .* k0.^2 - Nu.^2).*...
(I0(Nu) .* I0_2(Nu) + I2(Nu) .* I2_2(Nu)) - Nu.^2 .* 1./2 .* (I0(Nu) .* I2_2(Nu) + I0_2(Nu) .* I2(Nu))) .* Nu ;

del = 0.01*k0;
% del = eps;

wpoints = [(-1-1j).*del, (2+1j).*del k0+del+1j.*del k0+del-1j.*del del-1j.*del];

% wpoints = [(-1-1j).*del, (2+1j).*del];

% wpoints = [-del-1j*del, 2.*del+1j*del];

L = 1:100:3000;
for i = 1:length(L)
    E(i) =  integral(I, -del-1j.*L(i)*k0, del-1j.*L(i)*k0, 'Waypoints', wpoints) ./ sqrt(1 - (xmn./(k0.*R)).^2) .* 1./(120*pi);
%     Exx(i) = integral(Ixx, -L(i)*beta_rhop1, L(i)*beta_rhop1, 'Waypoints', wpoints);
end

figure;
plot(1:100:3000, real(E));
hold on;
plot(1:100:3000, imag(E));
grid on;



%% Function behavior
L = 1:100:2000;

for l = 1:length(L)
    
Nu = linspace(-L(l).*k0, L(l).*k0, 100*l);
% Nu = linspace(-eps-1j.*L(l).*k0, eps-1j.*L(l).*k0, 100*l);

dnu = Nu(2) - Nu(1);

kz = -1j .* sqrt(-(k0.^2 - (Nu).^2));

for i = 1:length(Nu)

        I0(i) = Lommel(0, R, beta_rhop, Nu(i), m - 1, m - 1);
        I2(i) = Lommel(0, R, beta_rhop, Nu(i), m + 1, m + 1);
        
        I0_2(i) = Lommel(0, R, beta_rhop1, Nu(i), m - 1, m - 1);
        I2_2(i) = Lommel(0, R, beta_rhop1, Nu(i), m + 1, m + 1);
end

Ixx = sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* (1/(4) .* omega .* mu) .* pi/2 .* (2.*k0^2 - Nu.^2)/kz .* I2 .* I2_2;

Exx(l) = sum(Ixx) .* dnu;


end

figure;
plot(L, abs(Exx), 'Linewidth', 2);

grid on;
% % % 
% % % 
% % % 
% % % plot(Nu./beta_rhop1, abs(Ixx), 'linewidth', 2); grid on;
