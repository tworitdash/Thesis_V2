mode = "TE";


[ficname,pathname] = uigetfile('*.efe','fichier ''.efe'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = load('result.txt');

rho_f = A(1:100, 1);
phi_f = A(1:100:36000, 2) * pi/180;
z_f = A(1, 1);
aux = A(:, 4:9);

[ficname,pathname] = uigetfile('*.hfe','fichier ''.hfe'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A_h = load('result.txt');

z_f_h = A_h(1, 1);
aux_h = A_h(:, 4:9);


[rho_f, phi_f] = meshgrid(rho_f, phi_f);

x = rho_f .* cos(phi_f);
y = rho_f .* sin(phi_f);

% f = 20*36000:1:21*36000 - 1;

f = 1:1:36000;

E_rho = aux(f, 1) + 1j .* aux(f, 2);
E_phi = aux(f, 3) + 1j .* aux(f, 4);
E_rho_reshape = reshape(E_rho, 100, 360);
E_phi_reshape = reshape(E_phi, 100, 360);

H_rho = aux_h(f, 1) + 1j .* aux_h(f, 2);
H_phi = aux_h(f, 3) + 1j .* aux_h(f, 4);
H_rho_reshape = reshape(H_rho, 100, 360);
H_phi_reshape = reshape(H_phi, 100, 360);


Ex = cos(phi_f) .* E_rho_reshape.' - sin(phi_f) .* E_phi_reshape.';
Ey = sin(phi_f) .* E_rho_reshape.' + cos(phi_f) .* E_phi_reshape.';


Hx = cos(phi_f) .* H_rho_reshape.' - sin(phi_f) .* H_phi_reshape.';
Hy = sin(phi_f) .* H_rho_reshape.' + cos(phi_f) .* H_phi_reshape.';


E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2);% + abs(E_z).^2);

E = reshape(E_tot, 100, 360);


H_tot = sqrt(abs(H_rho).^2 + abs(H_phi).^2);% + abs(E_z).^2);

H = reshape(H_tot, 100, 360);

e = 10;

if mode == "TM"

figure;
% h = pcolor(x,y, abs(H)./max(abs(H)));
h = pcolor(x,y, abs(H.')./max(max(abs(H.'))));
set(h,'ZData',-1+zeros(size(H.')))
hold on;

shading flat;

quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),real(Ex(1:e:end, 1:e:end)),real(Ey(1:e:end, 1:e:end)), 2, 'w');
colorbar;
colormap('jet');

xlabel('x [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y [m]', 'FontSize', 12, 'FontWeight', 'bold');

% title([mode, '_{',num2str(m),num2str(n),'}Linear Scale'], 'FontSize', 12, 'FontWeight', 'bold');
% print([mode, '_', num2str(m), num2str(n)], '-depsc');
title('Field data in Linear scale', 'FontSize', 12, 'FontWeight', 'bold');

% set(gcf,'renderer','Painters')
% print -depsc -tiff -r300 -painters TM11_Feko.eps;

else
    
figure;
h = pcolor(x,y, abs(E.')./max(max(abs(E.'))));
% surface(x,y, abs(E)); shading flat;
set(h,'ZData',-1+zeros(size(E.')))
hold on;
% 
shading flat;

quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),imag(Hx(1:e:end, 1:e:end)),imag(Hy(1:e:end, 1:e:end)), 2, 'w');
colorbar;
colormap('jet');

xlabel('x [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y [m]', 'FontSize', 12, 'FontWeight', 'bold');

% title([mode, '_{', num2str(m), num2str(n), '} Linear Scale'], 'FontSize', 12, 'FontWeight', 'bold');
title('Field data in Linear scale', 'FontSize', 12, 'FontWeight', 'bold');
% print([mode, '_', num2str(m), num2str(n)], '-depsc');

end

system('rm result.txt');