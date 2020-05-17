%% 

nomfic = '/Users/tworitdash/course/feko/NF_test_TE01.efe';

[E_rho_1, E_phi_1, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE11.efe';

[E_rho_2, E_phi_2, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE21.efe';

[E_rho_3, E_phi_3, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE31.efe';

[E_rho_4, E_phi_4, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE41.efe';

[E_rho_5, E_phi_5, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE12.efe';

[E_rho_6, E_phi_6, x_f, y_f] = FEKO_E(nomfic);

%% 
c0 = 3e8;

F = 14e9;
lamb = c0/F;

N = 6;

R = 2.3e-2;

z = 0.005; % + lamb/4;

er = 1;
mur = 1;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er .* er0;
mu = mur .* mu0;

drho = R/100;
dphi = pi/180;

[rho, phi] = meshgrid(eps:drho:R, eps:dphi:2*pi+eps);

%[Erho, Ephi, Ez] = E_r1(1:1:N, rho, phi, F, r, z, epsilon, mu);

[Erho, Ephi] = E_n([4 1 2 6 8 9], rho, phi, F, r, z, epsilon, mu, drho, dphi);

E_mat = [squeeze(max(max(abs(Erho(1, :, :))))) squeeze(max(max(abs(Erho(2, :, :))))) squeeze(max(max(abs(Erho(3, :, :))))) squeeze(max(max(abs(Erho(4, :, :))))) squeeze(max(max(abs(Erho(5, :, :))))) squeeze(max(max(abs(Erho(6, :, :))))); ...
    squeeze(max(max(abs(Ephi(1, :, :))))) (squeeze(max(max(abs(Ephi(2, :, :)))))) squeeze(max(max(abs(Ephi(3, :, :))))) squeeze(max(max(abs(Ephi(4, :, :))))) squeeze(max(max(abs(Ephi(5, :, :))))) squeeze(max(max(abs(Ephi(6, :, :)))))];

E_feko = [max(max(abs(E_rho_1))) max(max(abs(E_rho_2))) max(max(abs(E_rho_3))) max(max(abs(E_rho_4))) max(max(abs(E_rho_5))) max(max(abs(E_rho_6))); ...
    max(max(abs(E_phi_1))) max(max(abs(E_phi_2))) max(max(abs(E_phi_3))) max(max(abs(E_phi_4))) max(max(abs(E_phi_5))) max(max(abs(E_phi_6)))];

A = 1./(E_mat) .* E_feko;

[A_MM] = Norm([4 1 2 6 8 9]);

% figure;
% 
% surface(rho .* cos(phi), rho .* sin(phi), db(A(6) * abs(squeeze(Erho(6, :, :))))); shading flat;
% 
% figure;
% 
% surface(x_f, y_f, db(abs(E_rho_6).')); shading flat;

%% Investigation

% E_rho_feko = E_rho_1 + E_rho_2 + E_rho_3;
% E_phi_feko = E_phi_1 + E_phi_2 + E_phi_3;
% 
% E_rho_MM = squeeze(Erho(1, :, :)) .* A(1, 1) +  squeeze(Erho(2, :, :)) .* A(1, 2) +  squeeze(Erho(3, :, :)) .* A(1, 3);
% E_phi_MM = squeeze(Ephi(1, :, :)) .* A(2, 1) +  squeeze(Ephi(2, :, :)) .* A(2, 2) +  squeeze(Ephi(3, :, :)) .* A(2, 3);




m = [0 1 2 3 4];
n = [1 2];

[m, n] = meshgrid(m, n);

A_plot = [A(2, 1) A(1, 2) A(1, 3) A(1, 4) A(1, 5); 
          0 A(1, 6) 0 0 0];
figure;
surf(m, n, A_plot);


