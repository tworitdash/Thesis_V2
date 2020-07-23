[ficname,pathname] = uigetfile('*.dat','fichier ''.dat'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = readtable('result.txt');

th = A.Theta_deg_;
phi_ = eps:pi/180:2*pi+pi/180;
E_phi_0 = A.Total_Phi_0Deg__V_;
% E_phi_0 = A.Phi_0Deg_V_;
E_phi_90 = A.Total_Phi_90Deg__V_;
E_phi_45_th = A.Theta_Phi_45Deg__V_;
E_phi_45_ph = A.Phi_Phi_45Deg__V_;
% E_phi_90 = A.Phi_90Deg_V_;

figure(12);
hold on;
plot(th, db((abs(E_phi_0)/max(abs(E_phi_0)))), '-.', 'LineWidth', 2);
% plot(th, db((abs(E_phi_0./2.*k0./(4*pi)))), '-.', 'LineWidth', 2);
hold on;
plot(th, db((abs(E_phi_90)/max(abs(E_phi_90)))), '-.', 'LineWidth', 2);
% plot(th, db((abs(E_phi_90./2*k0/(4*pi)))), '-.', 'LineWidth', 2);


% xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('E_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
% grid on;
% % xlim([90 270]);
% ylim([-50 20]);

%% Co and Cross polar fields

Exp = 1./sqrt(2) .* (E_phi_45_th - E_phi_45_ph);
Eco = 1./sqrt(2) .* (E_phi_45_th + E_phi_45_ph);


figure(7);
hold on;
plot(th, db((abs(Eco)/max(abs(Eco)))), '-.', 'LineWidth', 2);
% plot(th, db((abs(E_phi_0./2.*k0./(4*pi)))), '-.', 'LineWidth', 2);
hold on;
plot(th, db((abs(Exp)/max(abs(Eco)))), '-.', 'LineWidth', 2);


figure(10);
hold on;
plot(th, db(abs(E_phi_45_th)./max(abs(E_phi_45_th))));
hold on;
plot(th, db(abs(E_phi_45_ph)./max(abs(E_phi_45_ph))));
hold on;

figure(10);
hold on;
plot(th, db(abs(E_phi_45_th)));
hold on;
plot(th, db(abs(E_phi_45_ph)));
hold on;