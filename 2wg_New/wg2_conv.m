clear;

F = 4e9:0.5e9:21e9;

rp = 0.0405319403216/2;   % radius of the bigger waveguide
rr = 0.0405319403216/2.1; % radius of the smaller waveguide

N = 30;

for i = 1:N

Nr = 1:1:i;
Np = 1:1:i;

erp = 1;
murp = 1;
err = 1;
murr = 1;


[X_til_pr] = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr);

for k = 1:length(F)
    disp('Freq iteration');
    disp(k);

[Spp_, Spr_, Srp_, Srr_] = GSM(Nr, Np, F(k), rp, rr, erp, murp, err, murr, X_til_pr);

    

    Spp(i).Spp_i(k, :, :) =  Spp_;
    Spr(i).Spr_i(k, :, :) =  Spr_;
    Srp(i).Srp_i(k, :, :) =  Srp_;
    Srr(i).Srr_i(k, :, :) =  Srr_;

end

% data5 = read(rfdata.data,'S_Feko_5modes_each.s10p');
% s_params_5 = extract(data5,'S_PARAMETERS');

figure(1);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Spp(i).Spp_i(:, 1, 1)))), 'LineWidth', 2); grid on;

figure(2);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Spr(i).Spr_i(:, 1, 1)))), 'LineWidth', 2); grid on;

figure(3);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(1, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Srp(i).Srp_i(:, 1, 1)))), 'LineWidth', 2); grid on;

figure(4);

% plot(F * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Srr(i).Srr_i(:, 1, 1)))), 'LineWidth', 2); grid on;

end
