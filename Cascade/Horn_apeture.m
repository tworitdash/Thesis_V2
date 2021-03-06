clear;


c0 = 3e8;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability


%% Feed dimensions

F = 5e9;
omega = 2 * pi * F;

lamb = c0./F;

k0 = 2.*pi./lamb;

rr = 2e-2; % Base radius
% rt = 4e-2; % Top redius
% rt = rr .* 3;
rt = 2 .* lamb;

n = 15; % number of transitions

% Length = 5e-2;
Length = 6 .* lamb;

R = linspace(rr, rt, n); % radius vector

drho = R(end)/100;
dphi = pi/180;

[rho, ph] = meshgrid(eps:drho:R(end), eps:dphi:2*pi+eps);

n = length(R);

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

epsilon = er .* er0;
mu = mur .* mu0;

L = ones(1, n) .* Length/n; % length of each waveguide section

% L(1) = 2.5 * L(1);
% L(end) = 1.25e-2;
L(1) = L(1) + lamb./4;

% [E_aperture_rho, E_aperture_phi] = Near_field_fun_2(er, mur, R, F, L, rho, ph, 20, drho, dphi);


[STT, STR, SRT, SRR, N] = GSM_N_SameAzimuth(R, L, er, mur, F, 20);

Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;
%% 
ModeNumberAper = N(end);


ap = zeros(N(end), 1);
ar = ones(N(1), 1);

if N(1) == 1
    STR_req = STR(1, 1:N(end));
else
    STR_req = squeeze(STR(1, 1:N(end), 1:N(1)));
end

bp =  STR_req * ar;  %+ squeeze(STT(1, 1:N(end), 1:N(end))) * ap;
Transmission_sum = bp;

%%

for i = 1:ModeNumberAper
    HigherModes = i+1:1:i+20;
    
    for d = 1:length(HigherModes)
    
        xmn(d) = str(HigherModes(d)).xmn;
        beta_rho(d) = xmn(d)./R(end);
        M(d) = str(HigherModes(d)).m;
        beta_z(d) = -1j .* sqrt(-(k0.^2 - (beta_rho(d)).^2));
        YTE(d) = beta_z(d)./(omega .* mu(end));
        ZTE(d) = 1./YTE(d);
    
    end
    
 for l = 1:length(HigherModes)
    for p = 1:length(HigherModes)
        disp('iteration inside');
        disp(p);
        disp('iteration outside');
        disp(l);
        
        if l == p
            Ymut(l, p) = Yin_Circular(HigherModes(l), HigherModes(p), k0, R(end), er(end), mur(end), 100) + YTE(p);
        else
            Ymut(l, p) = Yin_Circular(HigherModes(l), HigherModes(p), k0, R(end), er(end), mur(end), 100);
        end
    end
end

for e = 1:length(HigherModes)
    Y_rhs(e) = Yin_Circular(i, HigherModes(e), k0, R(end), er(end), mur(end), 100);
end

Dm(i, :) = Ymut\(-Y_rhs.');

Yii = Yin_Circular(i, i, k0, R(end), er(end), mur(end), 100);


xmnii = str(i).xmn;
beta_rhoii = xmnii./R(end);

beta_zii = -1j .* sqrt(-(k0.^2 - beta_rhoii.^2));
YTEii = beta_zii./(omega .* mu(end));
ZTEii = 1./YTEii;

yii(i) = Yii./YTEii;

yap(i) = (Yii + Dm(i, :) * Y_rhs.')./YTEii;

Gamma(i) = (1 - yap(i))./(1 + yap(i));

drho = R/100;
dph = pi/180;

[rho, ph] = meshgrid(eps:drho:R(end), eps:dph:2*pi);



NTEii_FEKO = 1j.* beta_zii ./ (beta_rhoii).^2 .* ZTEii;

ETE_i_rho(i, :, :) = NTEii_FEKO .* 1./rho .* besselj(1, beta_rhoii .* rho) .* sin(ph);
ETE_i_ph(i, :, :) = NTEii_FEKO .* beta_rhoii .* besselj_der(1, beta_rhoii .* rho) .* cos(ph);

% ETE_i = sqrt(abs(ETE_i_rho).^2 + abs(ETE_i_ph).^2);


HTE_i_ph(i, :, :) = squeeze(ETE_i_rho(i, :, :))./ZTEii;
HTE_i_rho(i, :, :) = -squeeze(ETE_i_ph(i, :, :))./ZTEii;

% HTE_fundamental = sqrt(abs(HTE_i_rho).^2 + abs(HTE_i_ph).^2);


ETE_higher_rho = zeros(size(rho));
ETE_higher_ph = zeros(size(rho));

HTE_higher_rho = zeros(size(rho));
HTE_higher_ph = zeros(size(rho));

for b = 1:length(HigherModes)

    NTE_FEKO(b) = 1j.* beta_z(b) ./ (beta_rho(b)).^2 .* ZTE(b);
%     NTE_FEKO(b) = sqrt((pi*(1)/2 .* (xmn(b).^2 - M(b).^2) .* (besselj(1, xmn(b))).^2).^(-1));
    ETE_higher_rho_b = (Dm(i, b)) .* NTE_FEKO(b) .* M(b)./rho .* besselj(M(b), beta_rho(b) .* rho) .* sin(M(b).*ph);
    ETE_higher_ph_b = (Dm(i, b)) .* NTE_FEKO(b) .* beta_rho(b) .* besselj_der(M(b), beta_rho(b) .* rho) .* cos(M(b).*ph);
    
    ETE_higher_rho = ETE_higher_rho + ETE_higher_rho_b;
    ETE_higher_ph = ETE_higher_ph + ETE_higher_ph_b;
    
    HTE_higher_ph = HTE_higher_ph - ETE_higher_rho_b./ZTE(b);
    HTE_higher_rho = HTE_higher_rho + ETE_higher_ph_b./ZTE(b);
    
end

ETE_ap_rho(i, :, :) = (1 + Gamma(i)) .* (squeeze(ETE_i_rho(i, :, :)));% + ETE_higher_rho);
ETE_ap_ph(i, :, :) = (1 + Gamma(i)) .* (squeeze(ETE_i_ph(i, :, :)));% + ETE_higher_ph);

HTE_ap_rho(i, :, :) = (1 - Gamma(i)) .* (squeeze(HTE_i_rho(i, :, :)));% + (1 + Gamma(i)) .* HTE_higher_rho;
HTE_ap_ph(i, :, :) = (1 - Gamma(i)) .* (squeeze(HTE_i_ph(i, :, :)));% + (1 + Gamma(i)) .* HTE_higher_ph;

end


E_inf_rho = zeros(size(rho));
E_inf_ph = zeros(size(rho));
E_fin_rho = zeros(size(rho));
E_fin_ph = zeros(size(rho));
H_inf_rho = zeros(size(rho));
H_inf_ph = zeros(size(rho));
H_fin_rho = zeros(size(rho));
H_fin_ph = zeros(size(rho));

%% 


for g = 1:ModeNumberAper
    
    E_inf_rho = E_inf_rho + squeeze(ETE_i_rho(g, :, :)) .* Transmission_sum(g);
    E_inf_ph = E_inf_ph + squeeze(ETE_i_ph(g, :, :)) .* Transmission_sum(g);
    
    H_inf_rho = H_inf_rho - squeeze(HTE_i_rho(g, :, :)) .* Transmission_sum(g);
    H_inf_ph = H_inf_ph - squeeze(HTE_i_ph(g, :, :)) .* Transmission_sum(g);
    
    E_fin_rho = E_fin_rho + squeeze(ETE_ap_rho(g, :, :)) .* Transmission_sum(g);
    E_fin_ph = E_fin_ph + squeeze(ETE_ap_ph(g, :, :)) .* Transmission_sum(g);
    
    H_fin_rho = H_fin_rho + squeeze(HTE_ap_rho(g, :, :)) .* Transmission_sum(g);
    H_fin_ph = H_fin_ph + squeeze(HTE_ap_ph(g, :, :)) .* Transmission_sum(g);
end


E_inf_ap = sqrt(abs(E_inf_rho).^2 + abs(E_inf_ph).^2);
E_fin_ap = sqrt(abs(E_fin_rho).^2 + abs(E_fin_ph).^2);

H_inf_ap = sqrt(abs(H_inf_rho).^2 + abs(H_inf_ph).^2);
H_fin_ap = sqrt(abs(H_fin_rho).^2 + abs(H_fin_ph).^2);



x = rho.*cos(ph);
y = rho.*sin(ph);


figure;
surface(x, y, db(abs(E_fin_ap))); shading flat;
colorbar;
colormap('jet');
xlabel('x[m]');
xlabel('y[m]');
title('E field on the aperture of a finite length circular waveguide');

figure;
surface(x, y, db(abs(H_fin_ap))); shading flat;
colorbar;
colormap('jet');
xlabel('x[m]');
ylabel('y[m]');
title('H field on the aperture of a finite length circular waveguide');

figure;
surface(x, y, db(abs(E_inf_ap))); shading flat;
xlabel('x[m]');
ylabel('y[m]');
title('E field on the aperture of an infinitely long circular waveguide');
colorbar;
colormap('jet');

figure;
surface(x, y, db(abs(H_inf_ap))); shading flat;
xlabel('x[m]');
ylabel('y[m]');
title('H field on the aperture of an infinitely long circular waveguide');
colorbar;
colormap('jet');


figure;
surface(x, y, ((E_inf_ap - E_fin_ap))); shading flat;
xlabel('x[m]');
ylabel('y[m]');
title('diff E Inf - fin');
colorbar;
colormap('jet');

figure;
surface(x, y, ((H_inf_ap - H_fin_ap))); shading flat;
xlabel('x[m]');
ylabel('y[m]');
title('diff H Inf - fin');
colorbar;
colormap('jet');

%% Far fields

dtheta = pi/180;
dphi = pi/180;

 
[theta, phi] = meshgrid(-pi/2+eps:dtheta:pi/2+eps, eps:dphi:2*pi+eps);

% [theta, phi] = meshgrid(eps:dtheta:pi+eps, eps:dphi:2*pi+eps);

Eth = zeros(size(theta));
Eph = zeros(size(theta));

for o = 1:ModeNumberAper
    
    HigherModes = o+1:1:o+20;
    
    [Eth_o, Eph_o] = FF_apertureFSCir(o, length(HigherModes), [1, Dm(o, :)], Gamma(o), theta, phi, F, er(end), mur(end), R(end));
    Eth = Eth + Eth_o .* Transmission_sum(o);
    Eph = Eph + Eph_o .* Transmission_sum(o);

end

% [Eth, Eph] = FF_apertureFSCir(1, 1, 1, Gamma(1), theta, phi, F, er(end), mur(end), R(end));
% Eth = Eth .* Transmission_sum(1);
% Eph = Eph + Eph .* Transmission_sum(1);


E_FF = sqrt(abs(Eth).^2 + abs(Eph).^2);

%% plots
figure(1);
hold on;
plot(theta(1, :).*180/pi, db(abs(E_FF(1, :))./max(abs(E_FF(1, :)))), 'LineWidth', 2);
hold on;
plot(theta(91, :).*180/pi, db(E_FF(91, :)./max(abs(E_FF(91, :)))), 'LineWidth', 2);

grid on;
% ylim([-50 30]);

%% Directivity
% theta = th;
dtheta = pi/180;
zeta = 120 .* pi;

U = abs(E_FF).^2./(2 .* zeta);

% P_rad_i = U(:, 91:end) .* sin(theta(:, 91:end)) .* dtheta .* dphi;

P_rad_i = U(:, 91:end) .* sin(theta(:, 91:end)) .* dtheta .* dphi;

P_rad = sum(sum(P_rad_i));

D = 4 .* pi .* U ./ P_rad;


figure;
plot(theta(1, :)*180/pi, db(D(1, :))/2, 'LineWidth', 2);
hold on;
plot(theta(91, :)*180/pi, db(D(91, :))/2, 'LineWidth', 2);

grid on;

%% 
Exp = (Eth(46, :) .* cos(theta(1, :)) - Eph(46, :));
Eco = (Eth(46, :) .* cos(theta(1, :)) + Eph(46, :));

figure(5);
hold on;
plot(theta(1, :).*180/pi, db(abs(Eco)./max(abs(Eco))), 'LineWidth', 2);
hold on;
plot(theta(91, :).*180/pi, db(abs(Exp)./max(abs(Eco))), 'LineWidth', 2);

