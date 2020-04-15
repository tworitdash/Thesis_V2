
nomfic = '/Users/tworitdash/course/feko/NF_test_TE11.efe';

[E_rho_1, E_phi_1, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TM01.efe';

[E_rho_2, E_phi_2, x_f, y_f] = FEKO_E(nomfic);

nomfic = '/Users/tworitdash/course/feko/NF_test_TE21.efe';

[E_rho_3, E_phi_3, x_f, y_f] = FEKO_E(nomfic);


E_rho_f = E_rho_1.'./max(max(E_rho_1.')) + E_rho_2.'./max(max(E_rho_2.')) + E_rho_3.'./max(max(E_rho_3.'));

E_phi_f = E_phi_1.'./max(max(E_phi_1.')) + E_phi_2.'./max(max(E_phi_2.')) + E_phi_3.'./max(max(E_phi_3.'));

% 
% 
% E_rho_f = E_rho_1.' + E_rho_2.' + E_rho_3.';
% 
% E_phi_f = E_phi_1.' + E_phi_2.' + E_phi_3.';



 E_abs_full = sqrt(abs(E_rho_f).^2 + abs(E_phi_f).^2);

% figure;
% 
% surface(x_f, y_f, db(abs(E_abs_full)./max(max((abs(E_abs_full)))))); shading flat;
% 
% colormap('jet');

figure;

surface(x_f, y_f, db(abs(E_rho_f)./max(max(abs(E_rho_f))))); shading flat;

% max(max(abs(E_rho_f)))

colormap('jet');

figure;

surface(x_f, y_f, db(abs(E_phi_f)./max(max(abs(E_phi_f))))); shading flat;

colormap('jet');


figure;

surface(x_f, y_f, angle((E_rho_f))); shading flat;

colormap('jet');

figure;

surface(x_f, y_f, angle((E_phi_f))); shading flat;
colormap('jet');

