% [ficname,pathname] = uigetfile('*.txt','fichier ''.efe'' a convertir ?');
% nomfic = [pathname ficname];
% i0 = find(ficname=='.');
% 
% system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = readtable('~/course/feko/cst_2/EX_Big_cst_im.txt');
x_ = A.Var1;
x = x_(1:47, 1);

y_ = A.Var2;
y = y_(1:47:end, 1);


aux = A(:, 4:9);

phi = linspace(-pi, pi, 47);


f = 1:1:2209;

E_x = sqrt(aux(f, 1).Var4.^2 + aux(f, 2).Var5.^2);
E_y = sqrt(aux(f, 3).Var6.^2 + aux(f, 4).Var7.^2);
E_z = sqrt(aux(f, 5).Var8.^2 + aux(f, 6).Var9.^2);

E_rho = sqrt(abs(E_x).^2 + abs(E_y).^2);
E_phi = atan(E_y./E_x);



E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2);% + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 47, 47);
E_rho_reshape = reshape(E_rho, 47, 47);
E_phi_reshape = reshape(E_phi, 47, 47);
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