%% 

nomfic = '/Users/tworitdash/course/feko/NF_test_TE11.efe';

[E_rho_1, E_phi_1, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TM01.efe';

[E_rho_2, E_phi_2, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE21.efe';

[E_rho_3, E_phi_3, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE01.efe';

[E_rho_4, E_phi_4, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TM11.efe';

[E_rho_5, E_phi_5, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE31.efe';

[E_rho_6, E_phi_6, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TM21.efe';

[E_rho_7, E_phi_7, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE41.efe';

[E_rho_8, E_phi_8, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE12.efe';

[E_rho_9, E_phi_9, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TM02.efe';

[E_rho_10, E_phi_10, x_f, y_f] = FEKO_E(nomfic);

%% 
c0 = 3e8;

F = 14e9;
lamb = c0/F;

N = 10;

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

[Erho, Ephi] = E_n(1:1:N, rho, phi, F, R, z, epsilon, mu, drho, dphi);
E_mat = zeros(2, length(N));

for i = 1:1:N
    if i == 4
        j = eps;
       
    elseif (i == 2) || (i == 10)
        k = eps;
    else
        j = 0;
        k = 0;
    end
    E_mat(1, i) = squeeze(max(max(abs(Erho(i, :, :)))));
    E_mat(2, i) = squeeze(max(max(abs(Ephi(i, :, :)))));
end


E_feko = [max(max(abs(E_rho_1))) max(max(abs(E_rho_2))) max(max(abs(E_rho_3)))...
    max(max(abs(E_rho_4))) max(max(abs(E_rho_5))) max(max(abs(E_rho_6)))...
    max(max(abs(E_rho_7))) max(max(abs(E_rho_8))) max(max(abs(E_rho_9))) max(max(abs(E_rho_10)));...
    max(max(abs(E_phi_1))) max(max(abs(E_phi_2))) max(max(abs(E_phi_3)))...
    max(max(abs(E_phi_4))) max(max(abs(E_phi_5))) max(max(abs(E_phi_6)))...
    max(max(abs(E_phi_7))) max(max(abs(E_phi_8))) max(max(abs(E_phi_9))) max(max(abs(E_phi_10)))];

A = 1./(E_mat) .* E_feko;

[A_MM] = Norm(1:1:N);


 


% figure;
% 
% surface(rho .* cos(phi), rho .* sin(phi), db(A(1) * abs(squeeze(Erho(1, :, :))))); shading flat;
% 
% figure;
% 
% surface(x_f, y_f, db(abs(E_rho_1).')); shading flat;

%% Investigation

E_rho_feko = E_rho_1 + E_rho_2 + E_rho_3 +  E_rho_4 +  E_rho_5 +  E_rho_6 +  E_rho_7 +  E_rho_8 +  E_rho_9 +  E_rho_10;
E_phi_feko = E_phi_1 + E_phi_2 + E_phi_3 + E_phi_4 +  E_phi_5 +  E_phi_6 +  E_phi_7 +  E_phi_8 +  E_phi_9 +  E_phi_10;

E_rho_MM = squeeze(Erho(1, :, :)) .* A(1, 1) +  squeeze(Erho(2, :, :)) .* A(1, 2) + ...
    squeeze(Erho(3, :, :)) .* A(1, 3) + squeeze(Erho(4, :, :)) .* A(2, 4) + ...
    squeeze(Erho(5, :, :)) .* A(1, 5) + squeeze(Erho(6, :, :)) .* A(1, 6) + ...
    squeeze(Erho(7, :, :)) .* A(1, 7) + squeeze(Erho(8, :, :)) .* A(1, 8) + ...
    squeeze(Erho(9, :, :)) .* A(1, 9) + squeeze(Erho(10, :, :)) .* A(1, 10);
E_phi_MM = squeeze(Ephi(1, :, :)) .* A(1, 1) +  squeeze(Ephi(2, :, :)) .* A(1, 2) + ...
    squeeze(Ephi(3, :, :)) .* A(1, 3) + squeeze(Ephi(4, :, :)) .* A(2, 4) + ...
    squeeze(Ephi(5, :, :)) .* A(1, 5) + squeeze(Ephi(6, :, :)) .* A(1, 6) + ...
    squeeze(Ephi(7, :, :)) .* A(1, 7) + squeeze(Ephi(8, :, :)) .* A(1, 8) + ...
    squeeze(Ephi(9, :, :)) .* A(1, 9) + squeeze(Ephi(10, :, :)) .* A(1, 10);

figure;

surface(rho .* cos(phi), rho .* sin(phi), db(E_rho_MM)); shading flat;
colormap('jet');

figure;

surface(x_f, y_f, db(abs(E_rho_feko).')); shading flat;
colormap('jet');

figure;

surface(rho .* cos(phi), rho .* sin(phi), db(E_phi_MM)); shading flat;
colormap('jet');

figure;

surface(x_f, y_f, db(abs(E_phi_feko).')); shading flat;
colormap('jet');


% figure;
% plot(1:1:N, A_MM);
% grid on;




