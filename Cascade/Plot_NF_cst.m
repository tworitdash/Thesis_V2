% [ficname,pathname] = uigetfile('*.txt','fichier ''.efe'' a convertir ?');
% nomfic = [pathname ficname];
% i0 = find(ficname=='.');
% 
% system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = readtable('~/course/feko/3modes_NF_CST_2wg.txt');
x_ = A.Var1;
x = x_(1:200, 1);

y_ = A.Var2;
y = y_(1:200:end, 1);

r = 0.0405319403216/2;

aux = A(:, 4:9);

phi = linspace(eps, 2*pi, 200);
rho = linspace(eps, r, 200);

[rho, phi] = meshgrid(rho, phi);


f = 1:1:40000;

E_x = sqrt(aux(f, 1).Var4.^2 + aux(f, 2).Var5.^2);
E_y = sqrt(aux(f, 3).Var6.^2 + aux(f, 4).Var7.^2);
E_z = sqrt(aux(f, 5).Var8.^2 + aux(f, 6).Var9.^2);


E_tot = sqrt(abs(E_x).^2 + abs(E_y).^2);% + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 200, 200);
E_x_reshape = reshape(E_x, 200, 200);
E_y_reshape = reshape(E_y, 200, 200);

E_rho_reshape = cos(phi) .* E_x_reshape + sin(phi) .* E_y_reshape;
E_phi_reshape = -sin(phi) .* E_x_reshape + cos(phi) .* E_y_reshape;

% E_z_reshape = reshape(E_z, 47, 47);


figure;

% surface(x_f, y_f, db((abs(E_tot_reshape).'))); shading flat;

surface(x, y, db((abs(E_tot_reshape).')./max(abs(E_tot_reshape).'))); shading flat;

colormap('jet');

figure;

surface(x, y, db((abs(E_rho_reshape).')./max(abs(E_rho_reshape).'))); shading flat;
colormap('jet');


figure;

surface(x, y, db((abs(E_phi_reshape).')./max(abs(E_phi_reshape).'))); shading flat;
colormap('jet');

figure;

surface(x, y, db((abs(E_x_reshape).')./max(abs(E_x_reshape).'))); shading flat;
colormap('jet');


figure;

surface(x, y, db((abs(E_y_reshape).')./max(abs(E_y_reshape).'))); shading flat;
colormap('jet');
