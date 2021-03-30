c0 = 3e8;

% R = [2e-2 4e-2];
R = [1.8e-2 3*1.8e-2];
% R = [1.8e-2 1.8e-2 3*1.8e-2 3*1.8e-2];
% Len = 5e-2;
Len = c0/5e9;

F = linspace(5e9, 21e9, 34);
k = 50;

% [STT, STR, SRT, SRR, num] = GSM_N_Vfinal(R, Len, F, 20);
[RLRR_TE11, SRR, SRT, STR, STT] = GSM_N_opt_allvar_V2_freq(R, Len, F, k);
figure;
plot(F .* 1e-9, db(abs(SRR(:, 1, 9))), 'LineWidth', 2); grid on;
figure;
plot(F .* 1e-9, db(abs(SRT(:, 1, 9))), 'LineWidth', 2); grid on;
figure;
hold on; plot(F .* 1e-9, db(abs(STR(:, 1, 1))), 'LineWidth', 2); grid on;
% figure;
% plot(F .* 1e-9, db(abs(STT(:, 1, 1))), 'LineWidth', 2); grid on;
% figure;
% plot(F .* 1e-9, angle(STR(:, 7, 7)), 'LineWidth', 2); grid on;


OPW.STT = STT;
OPW.STR = STR;
OPW.SRT = SRT;
OPW.SRR = SRR;
% fmin2.time_consumed = time_opt2;

save('Closed_Cone_lambda_length.mat', 'OPW');
%% FEKO results




data_feko = read(rfdata.data,'../../Thesis/Paper_Thesis/Cone/w_cone_lambda_length.s4p');
s_params_feko = extract(data_feko,'S_PARAMETERS');


figure;
plot(F .* 1e-9, (abs(SRT(:, 1, 9))), 'LineWidth', 2); grid on;

hold on;
% figure;
plot(F .* 1e-9, (abs(squeeze((s_params_feko(1, 4, :))))));
% hold on;

% plot(F .* 1e-9, (imag(squeeze((s_params_feko(1, 3, :))))));



%%

% m_feko = 1;
% n_feko = 3;
% m_mm = 1;
% n_mm = 1;
% 
% feko_str = squeeze(s_params_feko(m_feko, n_feko, :));
% mm_str = STR(:, m_mm, n_mm);
% 
% m_feko = 1;
% n_feko = 4;
% m_mm = 1;
% n_mm = 9;
% 
% feko_str_12 = squeeze(s_params_feko(m_feko, n_feko, :));
% mm_str_12 = STR(:, m_mm, n_mm);
% 
% m_feko = 1;
% n_feko = 1;
% m_mm = 1;
% n_mm = 1;
% 
% feko_srr = squeeze(s_params_feko(m_feko, n_feko, :));
% mm_srr = SRR(:, m_mm, n_mm);
% 
% m_feko = 1;
% n_feko = 2;
% m_mm = 1;
% n_mm = 9;
% 
% feko_srr_12 = squeeze(s_params_feko(m_feko, n_feko, :));
% mm_srr_12 = SRR(:, m_mm, n_mm);
% 
% 
% E_feko = sqrt(abs(feko_srr).^2 + abs(feko_str).^2 + abs(feko_str_12).^2 + abs(feko_srr_12).^2);
% 
% E_mm = sqrt(abs(mm_srr).^2 + abs(mm_str).^2 + abs(mm_str_12).^2 + abs(mm_srr_12).^2);
% 
% figure;
% plot(F .* 1e-9, E_mm, 'LineWidth', 2); grid on;
% 
% hold on;
% 
% plot(F .* 1e-9, E_feko, 'LineWidth', 2);
% 

%% Near Fields
data = load('Closed_Cone_lambda_length.mat');

STT = data.OPW.STT;
STR = data.OPW.STR;
SRT = data.OPW.SRT;
SRR = data.OPW.SRR;

R = [1.8e-2 3*1.8e-2];


F = linspace(5e9, 21e9, 34);

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

er = [1 1];
mur = [1 1];

epsilon = er .* er0;
mu = mur .* mu0;

% Number of modes propagating at eachport
% f_req = 16.15e9; % Frequency at which near fields are expected
f_req = 5e9;
n = length(R);
[rho, phi] = meshgrid(linspace(eps, R(end), 100), linspace(eps, 2*pi, 360));

drho = R/100;
dphi = pi/180;

z = c0/5e9/4 + c0/5e9 - 1e-3;

for i = 1:n
    f =  fc(R(i), er(i), mur(i));
    N_i  =  find(f < f_req);
    N(i) = length(N_i);
end

[Ep_rho, Ep_phi] = E_n(1:1:N(end), rho, phi, f_req, R(end), z, epsilon(end), mu(end), drho, dphi);


ap = zeros(N(end), 1);
ar = zeros(N(1), 1);
ar(1) = 1;

% if N(1) == 1
%     STR_req = STR(1, 1:9);
% else
    STR_req = squeeze(STR(1, 1:N(end), 1:N(1)));
% end
% STR_req(1, 1) = 0.4958 -1j .* 0.0802;

bp = squeeze(STT(1, 1:N(end), 1:N(end))) * ap + (STR_req * ar).';
Gamma_sum = ap + bp;

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
% E_aperture_z = zeros(size(rho));

for k = 1:1
    
    E_aperture_rho = E_aperture_rho + squeeze(Ep_rho(k, :, :)) .* (Gamma_sum(k)); %.* exp(1j .* pi/2);
    E_aperture_phi = E_aperture_phi + squeeze(Ep_phi(k, :, :)) .* (Gamma_sum(k)); %.* exp(1j .* pi/2);
%     E_aperture_z = E_aperture_z + squeeze(Ep_z(k, :, :)) .* (Gamma_sum(k));
end

E_res = sqrt(abs(E_aperture_rho).^2 + abs(E_aperture_phi).^2);

x = rho .* cos(phi);
y = rho .* sin(phi);
%% 
figure;

surface(x, y, (abs(E_aperture_rho))); shading flat;

colormap('jet');

figure;

surface(x, y, (abs(E_aperture_phi))); shading flat;

colormap('jet');

figure;

surface(x, y, (abs(E_res))); shading flat;

colormap('jet');
