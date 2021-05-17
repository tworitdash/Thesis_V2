clear;

c0 = 3e8;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

% l1 = c0/5e9/4;
% 
% L = [l1 l1];

Nr = 1:1:40;
Np = 1:1:40;
Nt = 1:1:41;

N = [Nr(end) Np(end) Nt(end)];

F = 4e9:0.5e9:21e9;

rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;
rt = 0.0405319403216/1.9;

R = [rr rp rt];
L = [0.001 0.02 0.001];

err = 1; murr = 1; erp = 1; murp = 1; ert = 1; murt = 1;

er = [err erp ert]; mur = [murr murp murt];

% F = linspace(5e9, 21e9, 35);

[X_til_rp] = Inner_prod(Nr, Np, rp, rr, erp, murp, err, murr);
[X_til_pt] = Inner_prod(Np, Nt, rt, rp, ert, murt, erp, murp);

parfor k = 1:length(F)

   [S11, S12, S21, S22] = GSM_V2(1:1:N(2), 1:1:N(3), F(k), R(3), R(2), er(3), mur(3), er(2), mur(2), X_til_pt);
   [S33, S34, S43, S44] = GSM_V2(1:1:N(1), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), X_til_rp);
   Sl = SL(R(2), F(k), 1:1:N(2), L(2));
   
   
   [STT_, STR_, SRT_, SRR_] = cascade_3(1:1:N(2), S11, S12, S21, S22, S33, S34, S43, S44, Sl);
     
   
   [slr] = SL(R(1), F(k), Nr, L(1));
   [slt] = SL(R(end), F(k), Nt, L(end));
   
   Stt(k, :, :) = slt * STT_ * slt; 
   Str(k, :, :) = slt * STR_ * slr; 
   Srt(k, :, :) = slr * SRT_ * slt; 
   Srr(k, :, :) = slr * SRR_ * slr;
end

save('Stt3_ratio_1_modes_20_new_analyt', 'Stt');
save('Str3_ratio_1_modes_20_new_analyt', 'Str');
save('Srt3_ratio_1_modes_20_new_analyt', 'Srt');
save('Srr3_ratio_1_modes_20_new_analyt', 'Srr');