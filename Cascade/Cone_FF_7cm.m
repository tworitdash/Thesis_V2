c0 = 3e8;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

%F = 21e9; % Frequency at which far field is requested
F = 7.5e9;
% fcutt = 20e9-1e-11;
% 
Str = load('Xmn.mat');
Xmn = Str.Xmn;

lamb = c0./F;

% rr = Xmn(1).xmn/(2 * pi * fcutt * sqrt(er0*mu0));
% rt = 2 * lamb;

rr = 2e-2; % Base radius
rt = 4e-2; % Top redius

n = 150; % number of transitions

Length = 7e-2;
% Length = 6 .* lamb;

R = linspace(rr, rt, n); % radius vector

% R = linspace(0.4583 .* lamb, 4.41 .* lamb, n);

drho = R(end)/100;
dphi = pi/180;

[rho, ph] = meshgrid(eps:drho:R(end), eps:dphi:2*pi+eps);

n = length(R);

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

epsilon = er .* er0;
mu = mur .* mu0;

L = ones(1, n) .* Length/n; % length of each waveguide section

% L(1) = lamb/4 + L(1);

L(1) = 2.5 * L(1);
L(end) = 1.25e-2;

% [E_aperture_rho, E_aperture_phi, E_aperture_z] = Near_field_fun(er, mur, R, F, L, rho, ph, 0);
% A = load('A.mat');
% A = A.A;

[E_aperture_rho, E_aperture_phi] = Near_field_fun_2(er, mur, R, F, L, rho, ph, 20, drho, dphi);

%%
% [Ex, Ey, Ez, Hx, Hy, Hz] = FF_M(E_aperture_rho, E_aperture_phi, rho, ph, theta, phi, F, drho, dphi);
% 
% H_abs = sqrt(Hx.^2 + Hy.^2 + Hz.^2);
% 
% E_abs = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
%%


x = rho .* cos(ph);
y = rho .* sin(ph);
% 
% figure;
% surface(x, y, db((abs(E_aperture))./max(max(abs(E_aperture))))); shading flat;
% colormap('jet');
% figure;

E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2);


figure(601)
surface(x,y, db(abs(E_aperture)./max(max(abs(E_aperture))))); shading flat;
colormap('jet');

figure(602);

surface(x,y, db(abs(E_aperture_rho)./max(max(abs(E_aperture_rho))))); shading flat;
% surface(x,y, db(abs(E_aperture_rho))); shading flat;
colormap('jet');

figure(603);

surface(x,y, db(abs(E_aperture_phi)./max(max(abs(E_aperture_phi))))); shading flat;
% surface(x,y, db(abs(E_aperture_phi))); shading flat;
colormap('jet');


figure(604);

surface(x, y, angle((E_aperture_rho))); shading flat;

colormap('jet');

figure(605);

surface(x, y, angle((E_aperture_phi))); shading flat;
colormap('jet');


%% FF
[theta, phi] = meshgrid(-pi/2-eps:pi/180:pi/2-eps, eps:pi/180:2*pi+eps);


% [theta, phi] = meshgrid(eps:pi/180:2*pi+eps, eps:pi/180:2*pi+eps);

[Eth, Eph] = FF_M_V2(E_aperture_rho, E_aperture_phi, rho, ph, theta, phi, F, drho, dphi);

E_abs = sqrt(abs(Eth).^2 + abs(Eph).^2);

% patternCustom(E_abs,theta,phi);

%FF Plots


figure(71);

hold on;  
plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))/max(abs(E_abs(1, :)))), 'LineWidth', 2);
% figure(42)
hold on;
plot(theta(91, :)*(180/pi), db(abs(E_abs(91, :))/max(abs(E_abs(91, :)))), 'LineWidth', 2);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;

ylim([-50 0]);


%%
% figure(22);
% 
% hold on;   
% plot(phi(:, 90)*(180/pi), db(abs(H_abs(:, 90))/max(abs(H_abs(:, 90)))), 'LineWidth', 2);
% hold on;
% plot(phi(:, end)*(180/pi), db(abs(H_abs(:, end))/max(abs(H_abs(:, end)))), 'LineWidth', 2);
% 
% 
% xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('H_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Far magnetic field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\theta = 0', '\theta = 90'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
% grid on;
% 
% ylim([-50 0]);
% 
% figure(23);
% 
% hold on;  
% plot(phi(:, 90)*(180/pi), db(abs(E_abs(:, 90))), 'LineWidth', 2);
% hold on;
% plot(phi(:, end)*(180/pi), db(abs(E_abs(:, end))), 'LineWidth', 2);
% 
% 
% xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('E_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\theta = 0', '\theta = 90'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
% grid on;
% 
% ylim([-50 0]);


%% 

