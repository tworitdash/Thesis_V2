clear;

% F = 4e9:0.5e9:21e9;

% rp = 0.0405319403216/2;   % radius of the bigger waveguide
% rr = 0.0405319403216/2.1; % radius of the smaller waveguide

rp = 0.04;
rr = 0.02;

N = 30;

parfor i = 1:N

Nr = 1:1:i;
Np = 1:1:4 * i;

erp = 1;
murp = 1;
err = 1;
murr = 1;

f =  fc(rr, 1, 1);
F = f(i) + eps;

[X_til_pr] = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr);

% for k = 1:length(F)
    disp('Mode iteration');
    disp(i);

[Spp_, Spr_, Srp_, Srr_] = GSM(Nr, Np, F, rp, rr, erp, murp, err, murr, X_til_pr);

    

    Spp(i) =  Spp_(1, 1);
    Spr(i) =  Spr_(1, 1);
    Srp(i) =  Srp_(1, 1);
    Srr(i) =  Srr_(1, 1);

end

save('Spp_conv_N_ratio_2', 'Spp');
save('Spr_conv_N_ratio_2', 'Spr');
save('Srp_conv_N_ratio_2', 'Srp');
save('Srr_conv_N_ratio_2', 'Srr');
% data5 = read(rfdata.data,'S_Feko_5modes_each.s10p');
% s_params_5 = extract(data5,'S_PARAMETERS');

figure(1);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(1:1:N, db(abs((Spp))), 'LineWidth', 2); grid on;

figure(2);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(1:1:N, db(abs((Spr))), 'LineWidth', 2); grid on;

figure(3);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(1:1:N, db(abs((Srp))), 'LineWidth', 2); grid on;

figure(4);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(1:1:N, db(abs((Srr))), 'LineWidth', 2); grid on;
