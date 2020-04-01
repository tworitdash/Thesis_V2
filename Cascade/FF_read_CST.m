% [ficname,pathname] = uigetfile('*.txt','fichier ''.txt'' a convertir ?');
% nomfic = [pathname ficname];
% i0 = find(ficname=='.');
% 
% system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

E0 = readtable('../../../feko/2wg_FF_CST_phi0.txt');
E90 = readtable('../../../feko/2wg_FF_CST.txt');

th = E0.Var1;
E_phi_0 = 10.^(E0.Var3./20);
E_phi_90 = 10.^(E90.Var3./20);


hold on;
plot(th, db(abs(E_phi_0)/max(abs(E_phi_0))), 'LineWidth', 2);
hold on;
plot(th, db(abs(E_phi_90)/max(abs(E_phi_90))), 'LineWidth', 2);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;
xlim([90 270]);
ylim([-50 0]);