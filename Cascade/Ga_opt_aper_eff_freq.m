function [e_ap] = Ga_opt_aper_eff_freq(R_cone, Len, f, k, focal_length, d)

c0 = 3e8;
lamb = c0./3e9;

for l = 1:length(f)
    F = f(k);

% R1 = R_vec(1);
% Rend = R_vec(2);

% R = linspace(R1, Rend, round(Len/(lamb/10)));

% R = linspace(R1, Rend, E);
lamb_opt_freq = c0./F(end);

num = round(Len./(lamb_opt_freq./10));

n_R = length(R_cone);

N_axis = round(num./(n_R-1));

for p = 1:n_R-1
    R_(:, p) = linspace(R_cone(p), R_cone(p+1), N_axis);
end

R = reshape(R_, 1, size(R_, 1) .* size(R_, 2));

n = length(R);

l1 = lamb/4;
L = [l1 ones(1, n - 1) * Len/n];

zeta = 120 .* pi;


er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(eps:dth:pi/2+eps, eps:dph:2*pi+eps);

[ModeNumberAper, Transmission_sum]  = Feed_Gamma_Dm_opt_V2L(R, L, F, k, er, mur);

[Eth, Eph, ~, ~, CO, XP] = Feed_FF_Superposition_V2L(ModeNumberAper, th, ph, F, er, mur, R, Transmission_sum);


% [Eth, Eph]  = Feed_FF(rr, rt, n, F, th, ph);

% E_FF = sqrt(abs(Eth).^2 + abs(Eph).^2);

C_spillover = 1/(2 .* zeta);

f_pattern_square = abs(Eth).^2 + abs(Eph).^2;

U_feed = C_spillover .* f_pattern_square;

CO_XP_square = abs(CO).^2 + abs(XP).^2;
CO_XP_half = abs(CO).^2 + 1./2 .* abs(XP).^2;


%% 

for i = 1:length(d)
    
    drho = d(i)/200; dphi = pi/ 180;

    [rho, phi] = meshgrid(eps:drho:d(i)/2, eps:dphi:2*pi);
    
    theta_ = 2 * atan(rho/(2 * focal_length));
    
    f_hash = focal_length./d(i);
    
    % Spillover efficiency (eta_s)
    
    theta_0(i) = 2 * acot(4 * f_hash);
    
    % numerator
    Int_u_n = U_feed(:,th(1,:)<=theta_0(i)) .* sin(th(:,th(1,:)<=theta_0(i))) .* dth .* dph;
    n_f(i) = sum(sum(Int_u_n));
    
    %denominator
    Int_u_d = U_feed .* sin(th) .* dth .* dph;
    d_f(i) = sum(sum(Int_u_d));
        
        
    eta_s(i) = n_f(i)/d_f(i);
    
    % Taper efficiency (eta_t)
%     
%     Eth_ = Eth(:, th(1,:)<=theta_);
%     Eph_ = Eph(:, th(1,:)<=theta_);
    
    %[Eth_, Eph_, Eco_, Exp_, CO_, XP_] = Feed_FF_Superposition(ModeNumberAper, Gamma, Dm, theta_, phi, F, er, mur, R, Transmission_sum, 20);
    [Eth_, Eph_, Eco_, ~, ~, ~] = Feed_FF_Superposition_V2L(ModeNumberAper, theta_, phi, F, er, mur, R, Transmission_sum);
    
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
    
    
    Area(i) = pi .* (d(i)^2)/4;
    
    e_tp(i) =  (1/Area(i)) * etp_n(i)./etp_d(i);
    

%Illumination efficiency

C_illumination(i) = 2 .* (cot(theta_0(i)./2)).^2;

eta_ill_n_i = abs(CO(:,th(1,:)<=theta_0(i))) .* tan(th(:,th(1,:)<=theta_0(i))./2) .* dth .* dph;
eta_ill_n(i) = (sum(sum(eta_ill_n_i))).^2;
eta_ill_d_i = (CO_XP_half(:,th(1,:)<=theta_0(i))) .* sin(th(:,th(1,:)<=theta_0(i))) .* dth .* dph;
eta_ill_d(i) = sum(sum(eta_ill_d_i ));

eta_ill(i) = C_illumination(i) .* eta_ill_n(i)./eta_ill_d(i);

%Polarization efficiency

eta_pol_n(i) = sum(sum((abs(CO_XP_half(:,th(1,:)<=theta_0(i)))) .* sin(th(:,th(1,:)<=theta_0(i))) .* dth .* dph));
eta_pol_d(i) = sum(sum(CO_XP_square(:,th(1,:)<=theta_0(i)) .* sin(th(:,th(1,:)<=theta_0(i))) .* dth .* dph));

eta_pol(i) = eta_pol_n(i)./eta_pol_d(i);

end

e_ap(l) = -eta_s .* eta_pol .* e_tp;

end
end



