clear;


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability


rt = 0.0405319403216/1.9;
rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;
rd = 2.2e-2;
re = 2.3e-2;




R = [rr rp rt rd re]; % radius vector

n = length(R);

F = 14e9; % Frequency of operation

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

epsilon = er .* er0;
mu = mur .* mu0;

L = 1e-3 * [0.5 1 1 20 0.5]; % length of each waveguide section

[STT, STR, SRT, SRR, N] = GSM_N(R, L, er, mur, F);

%% 
[rho, phi] = meshgrid(eps:R(1)/100:R(1),  eps:pi/180:2*pi-eps);

z = 0;

[Er_rho, Er_phi, Er_z] = E_r(1:1:N(1), rho, phi, F, R(1), z, epsilon(1), mu(1));

ap = zeros(N(end), 1);
ar = ones(N(1), 1);


br = squeeze(SRT(1, :, :)) * ap + squeeze(SRR(1, :, :)) * ar;

Gamma_sum = ar + br;

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));

for k = 1:N(1)
    E_aperture_rho = E_aperture_rho + squeeze(Er_rho(k, :, :)) .* abs(Gamma_sum(k));
    E_aperture_phi = E_aperture_phi + squeeze(Er_phi(k, :, :)) .* abs(Gamma_sum(k));
    E_aperture_z = E_aperture_z + squeeze(Er_z(k, :, :));% .* abs(Gamma_sum(k));
end

% E_aperture = sqrt(squeeze(abs(Er_rho(5, :, :))).^2 + squeeze(abs(Er_phi(5, :, :))).^2 + squeeze(abs(Er_z(5, :, :))).^2);
E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2 + abs(E_aperture_z).^2);
% E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2); % + abs(E_aperture_z).^2);


x = rho .* cos(phi);
y = rho .* sin(phi);

figure;

h = pcolor(x,y, db(abs(E_aperture)./max(abs(E_aperture))));

set(h,'ZData',-1+zeros(size(E_aperture_rho)))
hold on;

shading interp;

figure;

h = pcolor(x,y, db(abs(E_aperture_rho)./max(abs(E_aperture_rho))));

set(h,'ZData',-1+zeros(size(E_aperture_rho)))
hold on;

shading interp;

figure;

h = pcolor(x,y, db(abs(E_aperture_phi)./max(abs(E_aperture_phi))));

set(h,'ZData',-1+zeros(size(E_aperture_rho)))
hold on;

shading interp;

figure;

h = pcolor(x,y, db(abs(E_aperture_z)./max(abs(E_aperture_z))));

set(h,'ZData',-1+zeros(size(E_aperture_rho)))
hold on;

shading interp;




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

E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2 + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 100, 360);

E_rho_reshape = reshape(E_rho, 100, 360);
E_phi_reshape = reshape(E_phi, 100, 360);
E_z_reshape = reshape(E_z, 100, 360);


figure;

surface(x_f, y_f, db((abs(E_tot_reshape).')./max(abs(E_tot_reshape).'))); shading flat;

figure;

surface(x_f, y_f, db((abs(E_rho_reshape).')./max(abs(E_rho_reshape).'))); shading flat;

figure;

surface(x_f, y_f, db((abs(E_phi_reshape).')./max(abs(E_phi_reshape).'))); shading flat;

figure;

surface(x_f, y_f, db((abs(E_z_reshape).')./max(abs(E_z_reshape).'))); shading flat;


figure;

polar(phi(:, 1), (abs(E_aperture(:, 30)/max(abs(E_aperture(:, 30))))));
hold on;

polar(phi_f(:, 1), (abs(E_tot_reshape(30, :).'/max(abs(E_tot_reshape(30, :).')))));

grid on;


%%  ----------------------------------------------------------------------------------------------------------

[rho, phi] = meshgrid(eps:R(end)/100:R(end), eps:pi/180:2*pi+eps);

z = 0;

[Ep_rho, Ep_phi, Ep_z] = E_r(1:1:N(end), rho, phi, F, R(end), z, epsilon(end), mu(end));


ap = zeros(N(end), 1);
ar = ones(N(1), 1);

bp = squeeze(STT(1, :, :)) * ap + squeeze(STR(1, :, :)) * ar;
Gamma_sum = ap + bp;

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));

for k = 1:N(end)
    E_aperture_rho = E_aperture_rho + squeeze(Ep_rho(k, :, :)) .* (Gamma_sum(k));
    E_aperture_phi = E_aperture_phi + squeeze(Ep_phi(k, :, :)) .* (Gamma_sum(k));
    E_aperture_z = E_aperture_z + squeeze(Ep_z(k, :, :)) .* (Gamma_sum(k));
end

E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2 + abs(E_aperture_z).^2);

% E_aperture = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2);

x = rho .* cos(phi);
y = rho .* sin(phi);
% 
figure;
surface(x, y, db((abs(E_aperture))./max(abs(E_aperture)))); shading flat;

figure;

h = pcolor(x,y, db(abs(E_aperture_rho)./max(E_aperture_rho)));

set(h,'ZData',-1+zeros(size(E_aperture_rho)))
hold on;

shading interp;

figure;

h = pcolor(x,y, db(abs(E_aperture_phi)./max(E_aperture_phi)));

set(h,'ZData',-1+zeros(size(E_aperture_rho)))
hold on;

shading interp;

figure;

h = pcolor(x,y, db(abs(E_aperture_z)./max(E_aperture_z)));

set(h,'ZData',-1+zeros(size(E_aperture_rho)))
hold on;

shading interp;

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
E_z = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2 + abs(E_z).^2);

% E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2); % + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 100, 360);

E_rho_reshape = reshape(E_rho, 100, 360);
E_phi_reshape = reshape(E_phi, 100, 360);
E_z_reshape = reshape(E_z, 100, 360);


figure;

surface(x_f, y_f, db((abs(E_tot_reshape).')./max(abs(E_tot_reshape).'))); shading flat;

figure;

surface(x_f, y_f, db((abs(E_rho_reshape).')./max(abs(E_rho_reshape).'))); shading flat;

figure;

surface(x_f, y_f, db((abs(E_phi_reshape).')./max(abs(E_phi_reshape).'))); shading flat;

figure;

surface(x_f, y_f, db((abs(E_z_reshape).')./max(abs(E_z_reshape).'))); shading flat;


figure;

plot(rho(1, :), db(abs(E_aperture(90, :)/max(abs(E_aperture(90, :))))));
hold on;

plot(rho_f(1, :), db(abs(E_tot_reshape(:, 90)'/max(abs(E_tot_reshape(:, 90)))')));

grid on;



figure;

polar(phi(:, 1), (abs(E_aperture(:, 30)/max(abs(E_aperture(:, 30))))));
hold on;

polar(phi_f(:, 1), (abs(E_tot_reshape(30, :).'/max(abs(E_tot_reshape(30, :).')))));

grid on;

system('rm result.txt');
