[ficname,pathname] = uigetfile('*.dat','fichier ''.dat'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = readtable('result.txt');

th = A.Theta_deg_;
E_phi_0 = A.Phi_0Deg_V_;
E_phi_90 = A.Phi_90Deg_V_;


hold on;
plot(th, db(abs(E_phi_0)/max(abs(E_phi_0))), 'LineWidth', 2);
hold on;
plot(th, db(abs(E_phi_90)/max(abs(E_phi_90))), 'LineWidth', 2);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;
% xlim([90 270]);
ylim([-50 0]);