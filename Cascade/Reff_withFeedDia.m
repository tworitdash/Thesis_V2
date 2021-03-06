%% One waveguide with reflector efficiency
clear;
c0 = 3e8;


theta_0 = 53 .* pi/180;
d = 15;

R = linspace(0.1, 10, 100);
er = 1;
mur = 1;

focal_length = d/4 * cot(theta_0/2);

zeta = 120 .* pi;

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(eps:dth:pi/2+eps, eps:dph:2*pi);

f =  fc(R(1), er, mur);
F = f(1)+eps;

for i = 1:length(R)
    
[Eth, Eph, Eco, Exp, CO, XP] = FF_apertureFSCir(1, 1, 1, 0, th, ph, F, er, mur, R(i));

C_spillover = 1/(2 .* zeta);

f_pattern_square = abs(Eth).^2 + abs(Eph).^2;

U_feed = C_spillover .* f_pattern_square;

CO_XP_square = abs(CO).^2 + abs(XP).^2;
CO_XP_half = abs(CO).^2 + 1./2 .* abs(XP).^2;

    
    drho = d/200; dphi = pi/ 180;

    [rho, phi] = meshgrid(eps:drho:d/2, eps:dphi:2*pi);
    
    theta_ = 2 * atan(rho/(2 * focal_length));
    
    f_hash = focal_length/d;
    
    % Spillover efficiency (eta_s)
    
    % numerator
    Int_u_n = U_feed(:,th(1,:)<=theta_0) .* sin(th(:,th(1,:)<=theta_0)) .* dth .* dph;
    n_f(i) = sum(sum(Int_u_n));
    
    %denominator
    Int_u_d = U_feed .* sin(th) .* dth .* dph;
    d_f(i) = sum(sum(Int_u_d));
        
        
    eta_s(i) = n_f(i)/d_f(i);
    
    % Taper efficiency (eta_t)
%     
%     Eth_ = Eth(:, th(1,:)<=theta_);
%     Eph_ = Eph(:, th(1,:)<=theta_);
    
    [Eth_, Eph_, Eco_, Exp_, CO_, XP_] = FF_apertureFSCir(1, 1, 1, 0, theta_, phi, F, er, mur, R(i));

    C = ((4 .* focal_length)./(4 .* focal_length.^2 + rho.^2));
    
    Erho_ = - Eth_ .* C;
    Eph_ = - Eph_ .* C;
    Ecoa = - Eco_ .* C;
%    
%     Int_etp_n_rho = Erho_ .* rho .* drho .* dphi;
%     etp_n_rho(i) = sum(sum(Int_etp_n_rho));
%     
%     Int_etp_n_phi = Eph_ .* rho .* drho .* dphi;
%     etp_n_phi(i) = sum(sum(Int_etp_n_phi));
%     
%     etp_n(i) = abs(etp_n_rho(i)).^2 + abs(etp_n_phi(i)).^2;
    etp_n(i) = abs(sum(sum(Ecoa .* rho .* drho .* dphi))).^2;
    
    E_abs_rpz = abs(Erho_).^2 + abs(Eph_).^2;
    Int_etp_d = E_abs_rpz .* rho .* drho .* dphi;
    etp_d(i) = sum(sum(Int_etp_d));
    
    
    Area(i) = pi .* (d^2)/4;
    
    e_tp(i) =  (1/Area(i)) * etp_n(i)./etp_d(i);
    

%Illumination efficiency

C_illumination(i) = 2 .* (cot(theta_0./2)).^2;

eta_ill_n_i = abs(CO(:,th(1,:)<=theta_0)) .* tan(th(:,th(1,:)<=theta_0)./2) .* dth .* dph;
eta_ill_n(i) = (sum(sum(eta_ill_n_i))).^2;
eta_ill_d_i = (CO_XP_half(:,th(1,:)<=theta_0)) .* sin(th(:,th(1,:)<=theta_0)) .* dth .* dph;
eta_ill_d(i) = sum(sum(eta_ill_d_i ));

eta_ill(i) = C_illumination(i) .* eta_ill_n(i)./eta_ill_d(i);

%Polarization efficiency

eta_pol_n(i) = sum(sum((abs(CO_XP_half(:,th(1,:)<=theta_0))) .* sin(th(:,th(1,:)<=theta_0)) .* dth .* dph));
eta_pol_d(i) = sum(sum(CO_XP_square(:,th(1,:)<=theta_0) .* sin(th(:,th(1,:)<=theta_0)) .* dth .* dph));

eta_pol(i) = eta_pol_n(i)./eta_pol_d(i);

end

e_ap = eta_s .* e_tp;

figure;
plot(R, eta_s, 'LineWidth', 2);
hold on;
plot(R, eta_pol, 'LineWidth', 2);
hold on; 
plot(R, e_tp, 'LineWidth', 2);
hold on;
% plot(d, eta_ill, 'LineWidth', 2);
hold on;
plot(R, e_ap, 'LineWidth', 2);

grid on;


