c0 = 3e8;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

l1 = c0/5e9/4;

L = [l1 l1];

Nr = 1:1:50;
Np = 1:1:50;

N = [Nr(end) Np(end)];

F = 4e9:0.5e9:21e9;

rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;

R = [rr rp];

err = 1; murr = 1; erp = 1; murp = 1;

% F = linspace(5e9, 21e9, 35);

[X_til] = Inner_p2(Nr, Np, rp, rr, erp, murp, err, murr);
%%
for f = 1:length(F)
    [STT_, STR_, SRT_, SRR_] = GSM_V2(Nr, Np, F(f), rp, rr, erp, murp, err, murr, X_til);

    slr = SL(R(1), F(f), 1:1:N(1), L(1));
    slt = SL(R(end), F(f), 1:1:N(end), L(end));

    Spp(f, :, :) = slt * STT_ * slt; 
    Spr(f, :, :) = slt * STR_ * slr; 
    Srp(f, :, :) = slr * SRT_ * slt; 
    Srr(f, :, :) = slr * SRR_ * slr;
end

save('Spp2_ratio_1_modes_20_new', 'Spp');
save('Spr2_ratio_1_modes_20_new', 'Spr');
save('Srp2_ratio_1_modes_20_new', 'Srp');
save('Srr2_ratio_1_modes_20_new', 'Srr');

figure; plot(F.*1e-9, db(abs(squeeze(STT(:, 1, 1)))), 'LineWidth', 2); grid on;