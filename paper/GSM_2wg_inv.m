
%% GSM for one junction

clear; 

F = linspace(5e9, 25e9, 100);
% 
% rp = 0.0405319403216/2; % radius of the waveguide
% rr = 0.0405319403216/2.1;
rp = 1.8e-2*3; % radius of the waveguide
rr = 1.8e-2;
R = [rr rp];
n = length(R);
c0 = 3e8;

l = c0/5e9/4;

L = [l l];

err = 1;
murr = 1;
erp = 1;
murp = 1;

er = [err erp];
mu = [murr murp];


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;
epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;

for i = 1:n
    f =  fc(R(i), er(i), mu(i));
    N_i  =  find(f < F(end));
    N_(i) = length(N_i);
end

if (N_(end) < 50)
    %N = [20 round(20 * (R(end)/R(1))^2)];
    N = N_ + 50;
else
    N = N_;
end

[X_til] = Inner_p2(1:N(1), 1:N(2), R(2), R(1), er(2), mu(2), er(1), mu(1));
%%
for k = 1:length(F)
    disp('frequency iteration')
    disp(k)
    [STT_, STR_, SRT_, SRR_] = GSM_V2(1:N(1), 1:N(2), F(k), R(2), R(1), er(2), mu(2), er(1), mu(1), X_til);

    slr = SL(R(1), F(k), 1:1:N(1), L(1));
    slt = SL(R(end), F(k), 1:1:N(end), L(end));

    STT(k, :, :) = slt * STT_ * slt; 
    STR(k, :, :) = slt * STR_ * slr; 
    SRT(k, :, :) = slr * SRT_ * slt; 
    SRR(k, :, :) = slr * SRR_ * slr;
end
%%
for k = 1:length(F)
    disp('frequency iteration')
    disp(k)
    f_base =  fc(R(1), er(1), mu(1));
    Num_modes_prop  =  find(f_base < F(k));
    RLRR_TE11(k) = db(abs(sum(SRR(1, 1:Num_modes_prop(end)))).^2)./2;
end

figure;
plot(F * 1e-9, (2*(RLRR_TE11)), 'LineWidth', 2); grid on; hold on;
%%
data5_feko = read(rfdata.data,'Paper_2WG.s20p');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');


figure;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(STT(:, 1, 1)))), 'LineWidth', 2); grid on;
% hold on;
% figure;
plot(F * 1e-9, db(abs(squeeze(SRR(:, 1, 1)))), 'LineWidth', 2); grid on; hold on;
% plot(F_cst * 1e-9, db(abs(squeeze(s_params_5_cst(7, 7, :))))/2, 'LineWidth', 2); grid on;
% hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(1, 1, :)))), 'LineWidth', 2); grid on;


% xlim([9.5 21])

