
%% 
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

f = 1:1:36000;

E_z_1 = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

E_z_reshape_1 = reshape(E_z_1, 100, 360);

%% 

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


f = 1:1:36000;


E_z_2 = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

E_z_reshape_2 = reshape(E_z_2, 100, 360);

%% 

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


f = 1:1:36000;


E_z_3 = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

E_z_reshape_3 = reshape(E_z_3, 100, 360);


%% 
figure;

surface(x_f, y_f, db((abs(E_z_reshape_1).')./max(abs(E_z_reshape_1).'))); shading flat;
colormap('jet');


figure;

surface(x_f, y_f, db((abs(E_z_reshape_2).')./max(abs(E_z_reshape_2).'))); shading flat;
colormap('jet');

figure;

surface(x_f, y_f, db((abs(E_z_reshape_3).')./max(abs(E_z_reshape_3).'))); shading flat;
colormap('jet');

figure;

surface(x_f, y_f, db((abs(E_z_reshape_1 + E_z_reshape_2).')...
    ./max(abs(E_z_reshape_1 + E_z_reshape_2).'))); shading flat;
colormap('jet');
