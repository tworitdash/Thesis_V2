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
lamb = c0/F;

x = 0.6:0.01:1;

r = x .* lamb ./ 2;

for i = 1:length(r)

% R = 2e-2;
R = r(i);
% 
% for i = 1:length(f)
    
%     F = f(i);

beta = 2 * pi * F ./ c0;

omega = 2 * pi * F;




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
YTE = 1./ZTE;
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

I0 = @(Nu) Lommel2(0, R, beta_rhop, Nu, m - 1, m - 1);
I2 =  @(Nu) Lommel2(0, R, beta_rhop, Nu, m + 1, m + 1);
        
I0_2 = @(Nu) Lommel2(0, R, beta_rhop1, Nu, m - 1, m - 1);
I2_2 = @(Nu) Lommel2(0, R, beta_rhop1, Nu, m + 1, m + 1);
   
% end

% Ixx = @(Nu) beta_rhop .* beta_rhop1 .* (1/(4) .* omega .* mu) .* pi/2 .* (2.*k0^2 - Nu.^2)/kz(Nu) .* I2(Nu) .* I2_2(Nu) .* Nu;
% Iyy = @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* (1/(4) .* omega .* mu) .* pi/2 .* (2.*k0^2 - Nu.^2)/kz(Nu) .* I0(Nu) .* I0_2(Nu) .* Nu;

I =  @(Nu) (Nup) .* beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((2 .* k0.^2 - Nu.^2).*...
(I0(Nu) .* I0_2(Nu) + I2(Nu) .* I2_2(Nu)) - Nu.^2 .* 1./2 .* (I0(Nu) .* I2_2(Nu) + I0_2(Nu) .* I2(Nu))) .* Nu ;

del = 0.05.*k0;
% del = eps;

wpoints = [(-1-1j).*del, (2+1j).*del k0+del+1j.*del k0+del-1j.*del del-1j.*del];

% wpoints = [(-1-1j).*eps (2+1j).*eps beta_rhop-del+1j.*eps beta_rhop+1j.*del beta_rhop+del+1j.*eps k0-del+1j.*eps k0+1j.*del k0+del+1j.*eps k0+del - 1j.*del del-1j.*del];

wpoints2 = [(-1-1j).*del (2+1j).*del];

% wpoints = [-del-1j*del, 2.*del+1j*del];

L = 100;

% E =  integral(I, -del-1j.*L*k0, del-1j.*L*k0);%, 'Waypoints', wpoints); 
E_ = integral(I, eps, L.*k0);% 'Waypoints', wpoints2);
E_1 = integral(I, eps, L.*k0); % 'Waypoints', wpoints2);
E_2 = integral(I, eps, k0+eps);

Yin(i) = E_1./YTE;

Gamma_K(i) = (1 - Yin(i))./(1 + Yin(i));

end

figure;

plot(x, real(Yin), 'Linewidth', 2);
hold on;
plot(x, imag(Yin), 'Linewidth', 2);

grid on;

xlabel('2r/ \lambda', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Y_{in}', 'FontWeight', 'bold', 'FontSize', 16);
title('Input Admittance', 'FontWeight', 'bold', 'FontSize', 16);

legend({'g_{in}', 'b_{in}'}, 'location', 'northeast', 'FontWeight', 'bold', 'FontSize', 16);


figure(2);

plot(x, db(abs(Gamma_K)), 'Linewidth', 2);

grid on;

xlabel('2r/ \lambda', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('\Gamma in dB', 'FontWeight', 'bold', 'FontSize', 16);
title('Reflection coefficient', 'FontWeight', 'bold', 'FontSize', 16);


figure(3);

plot(x, (angle(Gamma_K))*180/pi, 'Linewidth', 2);

grid on;

xlabel('2r/ \lambda', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Phase \Gamma in deg', 'FontWeight', 'bold', 'FontSize', 16);
title('Reflection coefficient', 'FontWeight', 'bold', 'FontSize', 16);



% end

%% 
% 
% Str = load('Xmn.mat');
% 
% str = Str.Xmn;
% 
% xmn = str(9).xmn;
% 
% R = 2e-2;
% Omega = linspace(-100.*xmn./R, 100.*xmn./R, 201);
% % Omega = xmn./R;
% 
% for i = 1:length(Omega)
% 
%     I00(i) = Zeta(R, xmn/R, Omega(i), 0, 0);
%     I22(i) = Zeta(R, xmn/R, Omega(i), 2, 2);
% 
% end
% figure;
% plot(Omega/(xmn./R), abs(I00), 'LineWidth', 2);
% hold on;
% plot(Omega/(xmn./R), abs(I22), 'LineWidth', 2);
% grid on;
% 
% % figure;
% % plot(Omega/(xmn/R), abs(I22 + I00), 'LineWidth', 2);
% % hold on;
% % plot(Omega/(xmn/R), abs(I22 - I00), 'LineWidth', 2);
% 
% % 
% % grid on;
