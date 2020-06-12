clear;

c0 = 3e8; % speed of light
R = 2e-2; % Radius of the waveguide

er = 1;
mur = 1;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

F = 5e9;

lamb = c0/F;
omega = 2 * pi * F;

k0 = omega./c0;

L = 100;

N = 2:1:30;

Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;

for i = 1:length(N)
    
    xmn(i) = str(N(i)).xmn;
    beta_rho(i) = xmn(i)./R;
    M(i) = str(N(i)).m;
    beta_z(i) = -1j .* sqrt(-(k0.^2 - (beta_rho(i)).^2));
    YTE(i) = beta_z(i)./(omega .* mu);
    ZTE(i) = 1./YTE(i);
    
end

for l = 1:length(N)
    for p = 1:length(N)
        if l == p
            Ymut(l, p) = Yin_Circular(N(l), N(p), k0, R, er, mur, L) + YTE(p);
        else
            Ymut(l, p) = Yin_Circular(N(l), N(p), k0, R, er, mur, L);
        end
    end
end

for l = 1:length(N)
    Y_rhs(l) = Yin_Circular(1, N(l), k0, R, er, mur, L);
end

Dm = Ymut\(-Y_rhs.');

Y11 = Yin_Circular(1, 1, k0, R, er, mur, L);

xmn11 = str(1).xmn;
beta_rho11 = xmn11./R;

beta_z11 = -1j .* sqrt(-(k0.^2 - beta_rho11.^2));
YTE11 = beta_z11./(omega .* mu);
ZTE11 = 1./YTE11;

yap = (Y11 + Dm.' * Y_rhs.')./YTE11;

Gamma = (1 - yap)./(1 + yap);


%%


drho = R/100;
dph = pi/180;

[rho, ph] = meshgrid(eps:drho:R, eps:dph:2*pi);



NTE11_FEKO = 1j.* beta_z11 ./ (beta_rho11).^2 .* ZTE11;

ETE_fundamental_rho = NTE11_FEKO .* 1./rho .* besselj(1, beta_rho11 .* rho) .* sin(ph);
ETE_fundamental_ph = NTE11_FEKO .* beta_rho11 .* besselj_der(1, beta_rho11 .* rho) .* cos(ph);

ETE_fundamental = sqrt(abs(ETE_fundamental_rho).^2 + abs(ETE_fundamental_ph).^2);


HTE_fundamental_ph = ETE_fundamental_rho./ZTE11;
HTE_fundamental_rho = -ETE_fundamental_ph./ZTE11;

HTE_fundamental = sqrt(abs(HTE_fundamental_rho).^2 + abs(HTE_fundamental_ph).^2);


ETE_higher_rho = zeros(size(rho));
ETE_higher_ph = zeros(size(rho));

HTE_higher_rho = zeros(size(rho));
HTE_higher_ph = zeros(size(rho));

for b = 1:length(N)

    NTE_FEKO(b) = 1j.* beta_z(b) ./ (beta_rho(b)).^2 .* ZTE(b);
%     NTE_FEKO(b) = sqrt((pi*(1)/2 .* (xmn(b).^2 - M(b).^2) .* (besselj(1, xmn(b))).^2).^(-1));
    ETE_higher_rho_i = (Dm(b)) .* NTE_FEKO(b) .* M(b)./rho .* besselj(M(b), beta_rho(b) .* rho) .* sin(M(b).*ph);
    ETE_higher_ph_i = (Dm(b)) .* NTE_FEKO(b) .* beta_rho(b) .* besselj_der(M(b), beta_rho(b) .* rho) .* cos(M(b).*ph);
    
    ETE_higher_rho = ETE_higher_rho + ETE_higher_rho_i;
    ETE_higher_ph = ETE_higher_ph + ETE_higher_ph_i;
    
    HTE_higher_ph = HTE_higher_ph - ETE_higher_rho_i./ZTE(b);
    HTE_higher_rho = HTE_higher_rho + ETE_higher_ph_i./ZTE(b);
    
end



ETE_ap_rho = (1 + Gamma) .* (ETE_fundamental_rho + ETE_higher_rho);
ETE_ap_ph = (1 + Gamma) .* (ETE_fundamental_ph + ETE_higher_ph);

ETE_ap = sqrt(abs(ETE_ap_rho).^2 + abs(ETE_ap_ph).^2);

HTE_ap_rho = (1 - Gamma) .* (HTE_fundamental_rho) + (1 + Gamma) .* HTE_higher_rho;
HTE_ap_ph = (1 - Gamma) .* (HTE_fundamental_ph) + (1 + Gamma) .* HTE_higher_ph;

HTE_ap = sqrt(abs(HTE_ap_rho).^2 + abs(HTE_ap_ph).^2);

x = rho.*cos(ph);
y = rho.*sin(ph);




figure;
surface(x, y, db(abs(ETE_ap))); shading flat;
colorbar;
colormap('jet');
xlabel('x[m]');
xlabel('y[m]');
title('E field on the aperture of a finite length circular waveguide');

% figure;
% surface(x, y, db(abs(HTE_ap))); shading flat;
% colorbar;
% colormap('jet');
% xlabel('x[m]');
% ylabel('y[m]');
% title('H field on the aperture of a finite length circular waveguide');

figure;
surface(x, y, db(abs(ETE_fundamental))); shading flat;
xlabel('x[m]');
ylabel('y[m]');
title('E field on the aperture of an infinitely long circular waveguide');
colorbar;
colormap('jet');

% figure;
% surface(x, y, db(abs(HTE_fundamental))); shading flat;
% xlabel('x[m]');
% ylabel('y[m]');
% title('H field on the aperture of an infinitely long circular waveguide');
% colorbar;
% colormap('jet');

figure;
surface(x, y, (abs(ETE_fundamental) - abs(ETE_ap))); shading flat;
xlabel('x[m]');
ylabel('y[m]');
title('Diff E field (\infty long and finite length)');
colorbar;
colormap('jet');

% x = rho.*cos(ph);
% y = rho.*sin(ph);
% 
% figure;
% surface(x(1:360, :), y(1:360, :), db(abs(E_tot_reshape.') - abs(ETE_ap(1:360, :)))); shading flat;
% xlabel('x[m]');
% xlabel('y[m]');
% title('Diff E field FEKO and MATLAB');
% colorbar;
% colormap('jet');

%% Far fields

drho = R./100;
dphi = pi./180;
[theta, phi] = meshgrid(-pi/2-eps:pi/180:pi/2-eps, eps:pi/180:2*pi+eps);

[Eth, Eph] = FF_apertureFSCir(length(N)+1, [1, Dm.'], Gamma, theta, phi, F, er, mur, R);

E_abs = sqrt(abs(Eth).^2 + abs(Eph).^2);
% 
% [Eth_num, Eph_num] = FF_M_V2(ETE_ap_rho, ETE_ap_ph, rho, ph, theta, phi, F, drho, dphi);
% 
% E_abs_num = sqrt(abs(Eth_num).^2 + abs(Eph_num).^2);

% 
% [Eth_fund_ref, Eph_fund_ref] = FF_apertureFSCir(1, 1, Gamma, theta, phi, F, er, mur, R);
% 
% E_abs_fund_ref = sqrt(abs(Eth_fund_ref).^2 + abs(Eph_fund_ref).^2);
% 
[Eth_fund, Eph_fund] = FF_apertureFSCir(1, 1, 0, theta, phi, F, er, mur, R);

E_abs_fund = sqrt(abs(Eth_fund).^2 + abs(Eph_fund).^2);

% [Eth_fund_num, Eph_fund_num] = FF_M_V2(ETE_fundamental_rho, ETE_fundamental_ph, rho, ph, theta, phi, F, drho, dphi);
% 
% E_abs_fund_num = sqrt(abs(Eth_fund_num).^2 + abs(Eph_fund_num).^2);


figure(72);

hold on;  
plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))), '-.', 'LineWidth', 2);
hold on;
plot(theta(91, :)*(180/pi), db(abs(E_abs(91, :))), '-.', 'LineWidth', 2);

% 
% hold on;  
% plot(theta(1, :)*(180/pi), db(abs(E_abs_num(1, :))), '--', 'LineWidth', 2);
% hold on;
% plot(theta(91, :)*(180/pi), db(abs(E_abs_num(91, :))), '--', 'LineWidth', 2);

% hold on;  
% plot(theta(1, :)*(180/pi), db(abs(E_abs_fund_ref(1, :))), '--', 'LineWidth', 2);
% hold on;
% plot(theta(91, :)*(180/pi), db(abs(E_abs_fund_ref(91, :))), '--', 'LineWidth', 2);
% 
hold on;  
plot(theta(1, :)*(180/pi), db(abs(E_abs_fund(1, :))), '*', 'LineWidth', 1);
hold on;
plot(theta(91, :)*(180/pi), db(abs(E_abs_fund(91, :))), '*', 'LineWidth', 1);
hold on;  
% plot(theta(1, :)*(180/pi), db(abs(E_abs_fund_num(1, :))), '-.', 'LineWidth', 2);
% hold on;
% plot(theta(91, :)*(180/pi), db(abs(E_abs_fund_num(91, :))), '-.', 'LineWidth', 2);

xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\phi = 0 Finite length waveguide K space Rumsey', ...
    '\phi = 90 Finite length waveguide K space Rumsey', ...
    '\phi = 0 Infinite length waveguide analytical', ...
    '\phi = 90 Infinite length waveguide analytical', ...
    '\phi = 0 Finite length waveguide FEKO NF MATLAB FF', ...
    '\phi = 90 Finite length waveguide FEKO NF MATLAB FF', ...
    '\phi = 0 Infinite length waveguide FEKO NF MATLAB FF', ...
    '\phi = 90 Infinite length waveguide FEKO NF MATLAB FF'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;

ylim([-50 10]);


