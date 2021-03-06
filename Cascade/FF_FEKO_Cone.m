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

F = 6e9;

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

% 
figure;
surface(x_f, y_f, (abs(E_tot_reshape).')); shading flat; colorbar;
xlabel('x[m]');
ylabel('y[m]');
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
figure;

surface(x_f, y_f, abs((E_rho_reshape_).')); shading flat;

colormap('jet');

figure;

surface(x_f, y_f, abs((E_phi_reshape_).')); shading flat;
colormap('jet');
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
figure(37);

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


%% Directivity
dth = pi./180;
dph = pi/180;
zeta = 120 * pi;
U = abs(E_abs).^2 ./ (2 * zeta);
P_rad_i = U(1, 91:end) .* sin(theta(:, 91:end)) .* dth .* dph;
P_rad = sum(sum(P_rad_i));

D_ = 4*pi * U./P_rad;

figure(1);
hold on;
plot(theta(1, :)*180/pi, db(D_(1, :))/2, '*', 'LineWidth', 1);
grid on;
ylim([-40 25]);

% %% Aperture efficiency
% 
% dth = pi/180;
% dph = pi/180;
% 
% [th_, ph_] = meshgrid(eps:dth:pi/2+eps, eps:dph:2*pi+eps);
% 
% focal_length = 1; % Focal length of the reflector
% theta_0 = 50 .* pi./180; % Subtended angle given from SKA
% d = 4 .* focal_length .* tan(theta_0./2); % Reflectoor Diameter
% 
% 
% C_spillover = 1/(2 .* zeta);
% 
% f_pattern_square = abs(Eth).^2 + abs(Eph).^2;
% 
% U_feed = C_spillover .* f_pattern_square;
% 
% % CO_XP_square = abs(CO).^2 + abs(XP).^2;
% % CO_XP_half = abs(CO).^2 + 1./2 .* abs(XP).^2;
% 
% 
% %% 
% 
% for i = 1:length(d)
%     
% %     drho = d(i)/200; dphi = pi/ 180;
% % 
%     [rho, ph_o] = meshgrid(linspace(eps, R, 100), linspace(eps, 2*pi+eps, 360));
%     
%     theta_ = 2 * atan(rho/(2 * focal_length));
%     
%     
%     f_hash = focal_length./d(i);
%     
%     % Spillover efficiency (eta_s)
%     
%     theta_0(i) = 2 * acot(4 * f_hash);
%     
%     % numerator
%     Int_u_n = U_feed(:,th_(1,:)<=theta_0(i)) .* sin(th_(:,th_(1,:)<=theta_0(i))) .* dth .* dph;
%     n_f(i) = sum(sum(Int_u_n));
%     
%     %denominator
%     Int_u_d = U_feed(:, 91:end) .* sin(th_) .* dth .* dph;
%     d_f(i) = sum(sum(Int_u_d));
%         
%         
%     eta_s(i) = n_f(i)/d_f(i);
%     
%     % Taper efficiency (eta_t)
% %     
% %     Eth_ = Eth(:, th(1,:)<=theta_);
% %     Eph_ = Eph(:, th(1,:)<=theta_);
%     
%     %[Eth_, Eph_, Eco_, Exp_, CO_, XP_] = Feed_FF_Superposition(ModeNumberAper, Gamma, Dm, theta_, phi, F, er, mur, R, Transmission_sum, 20);
% %     [Eth_, Eph_, Eco_, ~, ~, ~] = Feed_FF_Superposition_V2L(ModeNumberAper, theta_, phi, F, er, mur, R, Transmission_sum);
%     [Eth_, Eph_] = FF_M_V2(E_rho_reshape_.', E_phi_reshape_.', rho, ph_o, theta_, ph_o, F, drho, dphi);
%     
%     Eco_ = sin(ph_o) .* Eth_ + cos(ph_o) .* Eph_;
% 
%     
%     C = ((4 .* focal_length)./(4 .* focal_length.^2 + rho.^2));
%     
%     Erho_ = - Eth_ .* C;
%     Eph_ = - Eph_ .* C;
%     Ecoa = - Eco_ .* C;
% %    
% %     Int_etp_n_rho = Erho_ .* rho .* drho .* dphi;
% %     etp_n_rho(i) = sum(sum(Int_etp_n_rho));
% %     
% %     Int_etp_n_phi = Eph_ .* rho .* drho .* dphi;
% %     etp_n_phi(i) = sum(sum(Int_etp_n_phi));
% %     
% %     etp_n(i) = abs(etp_n_rho(i)).^2 + abs(etp_n_phi(i)).^2;
%     etp_n(i) = abs(sum(sum(Ecoa .* rho .* drho .* dphi))).^2;
%     
%     E_abs_rpz = abs(Erho_).^2 + abs(Eph_).^2;
%     Int_etp_d = E_abs_rpz .* rho .* drho .* dphi;
%     etp_d(i) = sum(sum(Int_etp_d));
%     
%     
%     Area(i) = pi .* (d(i)^2)/4;
%     
%     e_tp(i) =  (1/Area(i)) * etp_n(i)./etp_d(i);
%     
% 
% %Illumination efficiency
% 
% % C_illumination(i) = 2 .* (cot(theta_0(i)./2)).^2;
% % 
% % eta_ill_n_i = abs(CO(:,th(1,:)<=theta_0(i))) .* tan(th(:,th(1,:)<=theta_0(i))./2) .* dth .* dph;
% % eta_ill_n(i) = (sum(sum(eta_ill_n_i))).^2;
% % eta_ill_d_i = (CO_XP_half(:,th(1,:)<=theta_0(i))) .* sin(th(:,th(1,:)<=theta_0(i))) .* dth .* dph;
% % eta_ill_d(i) = sum(sum(eta_ill_d_i ));
% % 
% % eta_ill(i) = C_illumination(i) .* eta_ill_n(i)./eta_ill_d(i);
% 
% %Polarization efficiency
% 
% eta_pol_n(i) = sum(sum((abs(Eco(:,th_(1,:)<=theta_0(i)))) .* sin(th_(:,th_(1,:)<=theta_0(i))) .* dth .* dph));
% eta_pol_d(i) = sum(sum(E_abs(:,th_(1,:)<=theta_0(i)) .* sin(th_(:,th_(1,:)<=theta_0(i))) .* dth .* dph));
% 
% eta_pol(i) = eta_pol_n(i)./eta_pol_d(i);
% 
% end
% 
% e_ap = eta_s .* eta_pol .* e_tp;
% 
% 
% 
% 
