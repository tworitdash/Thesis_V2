function [E_rho_reshape, E_phi_reshape, E_z_reshape, x_f, y_f] = E_FEKO_rec(nomfic)

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = load('result.txt');

x_f = A(1:1000, 1);
y_f = A(1:1000:1000000, 2);
z_f = A(1, 1);
aux = A(:, 4:9);

[x_f, y_f] = meshgrid(x_f, y_f);

% x_f = rho_f .* cos(phi_f);
% y_f = rho_f .* sin(phi_f);

% f = 20*36000:1:21*36000 - 1;

f = 1:1:1000000;

E_rho = sqrt(aux(f, 1).^2 + aux(f, 2).^2);
E_rho_ = aux(f, 1) + 1j .* aux(f, 2);
E_rho_reshape_ = reshape(E_rho_, 1000, 1000);

E_phi = sqrt(aux(f, 3).^2 + aux(f, 4).^2);
E_phi_ = aux(f, 3) + 1j .* aux(f, 4);
E_phi_reshape_ = reshape(E_phi_, 1000, 1000);

E_z = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

% E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2); % + abs(E_z).^2);

E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2);% + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 1000, 1000);
E_rho_reshape = reshape(E_rho, 1000, 1000);
E_phi_reshape = reshape(E_phi, 1000, 1000);
E_z_reshape = reshape(E_z, 1000, 1000);


% figure;
% 
% % surface(x_f, y_f, db((abs(E_tot_reshape).'))); shading flat;
% 
% surface(x_f, y_f, db((abs(E_tot_reshape).')./max(max(abs(E_tot_reshape).')))); shading flat;
% 
% colormap('jet');

% 
% % surface(x_f, y_f, (abs(E_tot_reshape).')); shading flat;
% 
% colormap('jet');

% figure;
% surface(x_f, y_f, db((abs(E_rho_reshape).'))); shading flat;
% % surface(x_f, y_f, db((abs(E_rho_reshape).')./max(max(abs(E_rho_reshape).')))); shading flat;
% colormap('jet');


% figure;
% surface(x_f, y_f, db((abs(E_phi_reshape).'))); shading flat;
% % surface(x_f, y_f, db((abs(E_phi_reshape).')./max(max(abs(E_phi_reshape).')))); shading flat;
% colormap('jet');
% colorbar;
% 
% figure;
% surface(x_f, y_f, db((abs(E_z_reshape).'))); shading flat;
% % surface(x_f, y_f, db((abs(E_z_reshape).')./max(max(abs(E_z_reshape).')))); shading flat;
% colormap('jet');

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

% 
% figure;
% 
% surface(x_f, y_f, angle((E_rho_reshape_).')); shading flat;
% 
% colormap('jet');
% 
% figure;
% 
% surface(x_f, y_f, angle((E_phi_reshape_).')); shading flat;
% colormap('jet');

system('rm result.txt');
end