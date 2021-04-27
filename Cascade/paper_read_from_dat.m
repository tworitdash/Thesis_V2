[ficname,pathname] = uigetfile('*.dat','fichier ''.dat'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = readtable('result.txt');
%
Theta = A.Theta_deg_;
E_matched = A.Phi_0Deg__;
E_open = A.Phi_90Deg__;

figure(3);
hold on;
plot(Theta, db(E_matched./max(E_matched)), '*');
hold on;
plot(Theta, db(E_open./max(E_matched)), '*');
%
% Theta = A.Theta_deg_;
% E_cof = A.LudwigIII_Co__Phi_0Deg__V_;
% E_xf = A.LudwigIII_Cross__Phi_45Deg__V_;
% 
% figure(81);
% hold on;
% plot(Theta, db(E_cof./max(E_cof)), '*');
% hold on;
% plot(Theta, db(E_xf./max(E_cof)), '*');
