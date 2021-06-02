function [E_rho_reshape_, E_phi_reshape_, x_f, y_f] = FEKO_E(nomfic)

system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  

A = load('result.txt');

rho_f = A(1:100, 1);
phi_f = A(1:100:36000, 2) * pi/180;
z_f = A(1, 1);
aux = A(:, 4:9);

[rho_f, phi_f] = meshgrid(rho_f, phi_f);

x_f = rho_f .* cos(phi_f);
y_f = rho_f .* sin(phi_f);

% f = 20*36000:1:21*36000 - 1;

f = 1:1:36000;

% E_rho = sqrt(aux(f, 1).^2 + aux(f, 2).^2);
E_rho_ = aux(f, 1) + 1j .* aux(f, 2);
E_rho_reshape_ = reshape(E_rho_, 100, 360);

% E_phi = sqrt(aux(f, 3).^2 + aux(f, 4).^2);
E_phi_ = aux(f, 3) + 1j .* aux(f, 4);
E_phi_reshape_ = reshape(E_phi_, 100, 360);



system('rm result.txt');
end