clear;
R = 1.93e-2;
F = 14e9;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

er = 1;
mur = 1;

epsilon = er .* er0;
mu = mur .* mu0;

[rho, phi] = meshgrid(eps:R/100:R,  eps:pi/180:2*pi-eps);

z = 0.0005;

N = 1:1:10;

[Er_rho, Er_phi, Er_z] = E_r(N, rho, phi, F, R, z, epsilon, mu);


E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));


for i = 1
    E_aperture_rho = E_aperture_rho + squeeze(Er_rho(i, :, :));
    E_aperture_phi = E_aperture_phi + squeeze(Er_phi(i, :, :));
    E_aperture_z = E_aperture_z + squeeze(Er_z(i, :, :));
end

E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2 + abs(E_aperture_z).^2);

% E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2);

x = rho .* cos(phi);
y = rho .* sin(phi);

figure;
surface(x, y, db((abs(E_aperture))./max(abs(E_aperture)))); shading flat;
% 
% figure;
% 
% h = pcolor(x,y, db(abs(E_aperture_rho)./max(E_aperture_rho)));
% 
% set(h,'ZData',-1+zeros(size(E_aperture_rho)))
% hold on;
% 
% shading interp;
% 
% figure;
% 
% h = pcolor(x,y, db(abs(E_aperture_phi)./max(E_aperture_phi)));
% 
% set(h,'ZData',-1+zeros(size(E_aperture_rho)))
% hold on;
% 
% shading interp;

figure;

surface(x,y, db(abs(E_aperture_z)./max(abs(E_aperture_z)))); shading flat;

% figure;
% 
% surface(x,y, (angle(E_aperture_z))); shading flat;
% 
% set(h,'ZData',-1+zeros(size(E_aperture_rho)))
% hold on;
% 
% shading interp;


[ficname,pathname] = uigetfile('*.efe','fichier ''.efe'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');
% <<<<<<< HEAD
% system(['"C:\Program Files (x86)\GnuWin32\bin\sed" -e "/^#/d;/^*/d" ',' "',nomfic,'"| "C:\Program Files (x86)\GnuWin32\bin\tr" -s " " " " > result.txt']);    
% =======
system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608
A = load('result.txt');

rho_f = A(1:100, 1);
phi_f = A(1:100:36000, 2) * pi/180;
z_f = A(1, 1);
aux = A(:, 4:9);

% <<<<<<< HEAD
[rho_f, phi_f] = meshgrid(rho_f, pi - phi_f);
% =======
% [rho_f, phi_f] = meshgrid(rho_f, phi_f);
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608

x_f = rho_f .* cos(phi_f);
y_f = rho_f .* sin(phi_f);

% f = 20*36000:1:21*36000 - 1;

f = 1:1:36000;

E_rho = sqrt(aux(f, 1).^2 + aux(f, 2).^2);
E_phi = sqrt(aux(f, 3).^2 + aux(f, 4).^2);

% E_z_vec = aux(f, 5) + 1j .* aux(f, 6);



%%
[ficname,pathname] = uigetfile('*.efe','fichier ''.efe'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');
% <<<<<<< HEAD
% system(['"C:\Program Files (x86)\GnuWin32\bin\sed" -e "/^#/d;/^*/d" ',' "',nomfic,'"| "C:\Program Files (x86)\GnuWin32\bin\tr" -s " " " " > result.txt']);    
% =======
system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608
A = load('result.txt');

rho_f = A(1:100, 1);
phi_f = A(1:100:36000, 2) * pi/180;
z_f = A(1, 1);
aux = A(:, 4:9);

%%

% <<<<<<< HEAD
[rho_f, phi_f] = meshgrid(rho_f, phi_f);
% =======
% [rho_f, phi_f] = meshgrid(rho_f, phi_f);
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608

x_f = rho_f .* cos(phi_f);
y_f = rho_f .* sin(phi_f);

E_z = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2 + abs(E_z).^2);

% E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2); % + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 100, 360);

E_rho_reshape = reshape(E_rho, 100, 360);
E_phi_reshape = reshape(E_phi, 100, 360);
E_z_reshape = reshape(E_z, 100, 360);
% E_z_vec_reshape = reshape(E_z_vec, 100, 360);


figure;

surface(x_f, y_f, db((abs(E_tot_reshape).')./max(abs(E_tot_reshape).'))); shading flat;

% figure;
% 
% surface(x_f, y_f, db((abs(E_rho_reshape).')./max(abs(E_rho_reshape).'))); shading flat;
% 
% figure;
% 
% surface(x_f, y_f, db((abs(E_phi_reshape).')./max(abs(E_phi_reshape).'))); shading flat;

figure;

surface(x_f, y_f, db((abs(E_z_reshape).')./max(abs(E_z_reshape).'))); shading flat;



% figure;
% 
% surface(x_f, y_f, (angle(E_z_vec_reshape).')); shading flat;


system('rm result.txt')