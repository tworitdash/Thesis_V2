

F = 4e9:0.5e9:21e9;

rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;

c0 = 3e8;

err = 1;
murr = 1;
erp = 1;
murp = 1;


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;
epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;


Nr = 1:1:15;
Np = 1:1:15;

X_til = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr);


parfor k = 1:length(F)


[Spp_, Spr_, Srp_, Srr_] = GSM(Nr, Np, F(k), rp, rr, erp, murp, err, murr, X_til);

Spp(k, :, :) = Spp_;
Spr(k, :, :) = Spr_;
Srp(k, :, :) = Srp_;
Srr(k, :, :) = Srr_;


end


F_cst = linspace(4e9, 21e9, 1001);
data5_cst = read(rfdata.data,'2wg_cst.s38p');
s_params_5_cst = extract(data5_cst,'S_PARAMETERS');


data5_feko = read(rfdata.data,'2wg_feko.s20p');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');

figure;

plot(F * 1e-9, db(abs(squeeze(Spp(:, 4, 4))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(25,25, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(4, 4, :))))/2, 'LineWidth', 2); grid on;


xlim([5 21])

