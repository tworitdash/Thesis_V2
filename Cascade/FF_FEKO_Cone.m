[ficname,pathname] = uigetfile('*.efe','fichier ''.efe'' a convertir ?');
% [ficname,pathname] = uigetfile('*.hfe','fichier ''.hfe'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = load('result.txt');

rho_f = A(1:100, 1);
phi_f = A(1:100:36000, 2) * pi/180;
z_f = A(1, 1);
aux = A(:, 4:9);

[rho_f, phi_f] = meshgrid(rho_f, phi_f);
R = rho_f(end);

x_f = rho_f .* cos(phi_f);
y_f = rho_f .* sin(phi_f);

F = 3.3250e9;

% f = 20*36000:1:21*36000 - 1;

f = 1:1:36000;

E_rho = sqrt(aux(f, 1).^2 + aux(f, 2).^2);
E_rho_ = aux(f, 1) + 1j .* aux(f, 2);
E_rho_reshape_ = reshape(E_rho_, 100, 360);

E_phi = sqrt(aux(f, 3).^2 + aux(f, 4).^2);
E_phi_ = aux(f, 3) + 1j .* aux(f, 4);
E_phi_reshape_ = reshape(E_phi_, 100, 360);

E_z = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

% E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2); % + abs(E_z).^2);

E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2);% + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 100, 360);
E_rho_reshape = reshape(E_rho, 100, 360);
E_phi_reshape = reshape(E_phi, 100, 360);
E_z_reshape = reshape(E_z, 100, 360);


% figure;
% 
% % surface(x_f, y_f, db((abs(E_tot_reshape).'))); shading flat;
% % 
% surface(x_f, y_f, db((abs(E_tot_reshape).')./max(max(abs(E_tot_reshape).')))); shading flat;
% % % 
% colormap('jet');

% % 
% figure;
% surface(x_f, y_f, db(abs(E_tot_reshape).')); shading flat; colorbar;
% xlabel('x[m]');
% ylabel('y[m]');
% 
% colormap('jet');
% 
% figure;
% 
% surface(x_f, y_f, db((abs(E_rho_reshape).')./max(max(abs(E_rho_reshape).')))); shading flat;
% colormap('jet');
% 
% 
% figure;
% 
% surface(x_f, y_f, db((abs(E_phi_reshape).')./max(max(abs(E_phi_reshape).')))); shading flat;
% colormap('jet');
% % 
% 
% figure;
% 
% surface(x_f, y_f, angle((E_rho_reshape_).')); shading flat;
% 
% colormap('jet');
% 
% figure;
% 
% surface(x_f, y_f, angle((E_phi_reshape_).')); shading flat;
% colormap('jet');
%%

drho = R/100;
dphi = pi/180;

[rho, ph] = meshgrid(linspace(eps, R, 100), linspace(eps, 2*pi+eps, 360));

[theta, phi] = meshgrid(-pi/2-eps:pi/180:pi/2-eps, eps:pi/180:2*pi+eps);
%
% [Ex, Ey, Ez, Hx, Hy, Hz] = FF_M(E_rho_reshape_.', E_phi_reshape_.', rho, ph, theta, phi, F, drho, dphi);
% 
% E_abs = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
% 
% H_abs = sqrt(Hx.^2 + Hy.^2 + Hz.^2);

%%

[Eth, Eph] = FF_M_V2(E_rho_reshape_.', E_phi_reshape_.', rho, ph, theta, phi, F, drho, dphi);
E_abs = sqrt(abs(Eth).^2 + abs(Eph).^2);

Exp = cos(phi) .* Eth - sin(phi) .* Eph;
Eco = sin(phi) .* Eth + cos(phi) .* Eph;

% patternCustom(E_abs,theta,phi);

% figure(20);
% % 
% hold on;   
% plot(theta(1, :)*(180/pi), db(abs(H_abs(1, :))/max(abs(H_abs(1, :)))), 'LineWidth', 2);
% hold on;
% plot(theta(90, :)*(180/pi), db(abs(H_abs(90, :))/max(abs(H_abs(90, :)))), 'LineWidth', 2);
% 
% 
% xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('H_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Far magnetic field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
% grid on;
% 
% ylim([-50 0]);
%%
figure(80);

hold on;
plot(theta(1, :)*(180/pi), db((abs(E_abs(1, :))/max(abs(E_abs(1, :))))), '*', 'LineWidth', 0.5);
% plot(theta(1, :)*(180/pi), db((abs(E_abs(1, :)))), 'o', 'LineWidth', 1);
% figure(42)
hold on;
plot(theta(90, :)*(180/pi), db((abs(E_abs(90, :))/max(abs(E_abs(90, :))))), '*', 'LineWidth', 0.5);
% plot(theta(90, :)*(180/pi), db((abs(E_abs(90, :)))), 'o', 'LineWidth', 1);

xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs} (dB) normalized ', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;

ylim([-40 5]);

figure(81);

hold on;
plot(theta(1, :)*(180/pi), db((abs(Eco(1, :))./max(max(abs(Eco))))), '*', 'LineWidth', 0.5);
% plot(theta(1, :)*(180/pi), db((abs(E_abs(1, :)))), 'o', 'LineWidth', 1);
% figure(42)
hold on;
plot(theta(90, :)*(180/pi), db((abs(Exp(45, :))./max(max(abs(Eco))))), '*', 'LineWidth', 0.5);
% plot(theta(90, :)*(180/pi), db((abs(E_abs(90, :)))), 'o', 'LineWidth', 1);

xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{co}, E_{xp} in (dB) normalized', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');
ylim([-40 5]);
grid on;
system('rm result.txt');