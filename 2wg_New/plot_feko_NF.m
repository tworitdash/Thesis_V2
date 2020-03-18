[ficname,pathname] = uigetfile('*.efe','fichier ''.efe'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = load('result.txt');

rho_f = A(1:100, 1);
phi_f = A(1:100:36000, 2) * pi/180;
z_f = A(1, 1);
aux = A(:, 4:9);

[rho_f, phi_f] = meshgrid(rho_f, pi - phi_f);

x_f = rho_f .* cos(phi_f);
y_f = rho_f .* sin(phi_f);

% f = 20*36000:1:21*36000 - 1;

f = 1:1:36000;

E_rho = sqrt(aux(f, 1).^2 + aux(f, 2).^2);
E_phi = sqrt(aux(f, 3).^2 + aux(f, 4).^2);
E_z = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

% E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2); % + abs(E_z).^2);

E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2);% + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 100, 360);
E_rho_reshape = reshape(E_rho, 100, 360);
E_phi_reshape = reshape(E_phi, 100, 360);
E_z_reshape = reshape(E_z, 100, 360);


% figure;

% surface(x_f, y_f, db((abs(E_tot_reshape).'))); shading flat;

% surface(x_f, y_f, db((abs(E_tot_reshape).')./max(abs(E_tot_reshape).'))); shading flat;
% 
% colormap('jet');

% 
% % surface(x_f, y_f, (abs(E_tot_reshape).')); shading flat;
% 
% colormap('jet');

% figure;
% 
% surface(x_f, y_f, db((abs(E_rho_reshape).')./max(abs(E_rho_reshape).'))); shading flat;
% colormap('jet');
% 
% 
% figure;
% 
% surface(x_f, y_f, db((abs(E_phi_reshape).')./max(abs(E_phi_reshape).'))); shading flat;
% colormap('jet');


figure;

surface(x_f, y_f, db((abs(E_z_reshape).')./max(abs(E_z_reshape).'))); shading flat;
colormap('jet');

%% cuts

% figure;
% 
% polar(phi(:, 1), (abs(E_aperture(:, 30)/max(abs(E_aperture(:, 30))))));
% hold on;
% 
% polar(phi_f(:, 1), (abs(E_tot_reshape(30, :).'/max(abs(E_tot_reshape(30, :).')))));
% 
% grid on;
% 
% figure;
% 
% polar(phi(:, 1), (abs(E_aperture(:, 30)/max(abs(E_aperture(:, 30))))));
% hold on;
% 
% polar(phi_f(:, 1), (abs(E_tot_reshape(30, :).'/max(abs(E_tot_reshape(30, :).')))));
% 
% grid on;


% figure;
% 
% polar(phi(:, 1), (abs(E_aperture(:, 30)/max(abs(E_aperture(:, 30))))));
% hold on;
% 
% polar(phi_f(:, 1), (abs(E_tot_reshape(30, :).'/max(abs(E_tot_reshape(30, :).')))));
% 
% grid on;

system('rm result.txt');