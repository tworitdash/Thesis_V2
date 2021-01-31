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

% F = 17e9;

% f = 20*36000:1:21*36000 - 1;

f = 1:1:36000;

E_rho_f = sqrt(aux(f, 1).^2 + aux(f, 2).^2);
E_rho_ = aux(f, 1) + 1j .* aux(f, 2);
E_rho_reshape_ = reshape(E_rho_, 100, 360);

E_phi_f = sqrt(aux(f, 3).^2 + aux(f, 4).^2);
E_phi_ = aux(f, 3) + 1j .* aux(f, 4);
E_phi_reshape_ = reshape(E_phi_, 100, 360);

% E_z = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

% E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2); % + abs(E_z).^2);

E_tot = sqrt(abs(E_rho_f).^2 + abs(E_phi_f).^2);% + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 100, 360);
E_rho_reshape = reshape(E_rho_f, 100, 360);
E_phi_reshape = reshape(E_phi_f, 100, 360);
% E_z_reshape = reshape(E_z, 100, 360);




