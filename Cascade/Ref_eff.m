clear;

c0 = 3e8;
zeta = 120 .* pi;
F = 5e9;
lamb = c0./F;

%% Reflector dimensions

focal_length = 1;
d = linspace(0.2, 2, 100);

%% Feed dimensions

rr = 2e-2;
rt = 1.5 .* lamb;
n = 15;
R = linspace(rr, rt, n);

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

Length = 0.5 .* lamb; % length of the conical waveguide

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(eps:dth:pi/2-eps, eps:dph:2*pi+eps);

[Gamma, Dm, ModeNumberAper, Transmission_sum]  = Feed_Gamma_Dm(R, Length, F, 5, er, mur);

[Eth, Eph] = Feed_FF_Superposition(ModeNumberAper, Gamma, Dm, th, ph, F, er, mur, R, Transmission_sum);


% [Eth, Eph]  = Feed_FF(rr, rt, n, F, th, ph);

E_FF = sqrt(abs(Eth).^2 + abs(Eph).^2);


%% 
figure; 
plot(th(1, :).*180/pi, db(E_FF(1, :)), 'LineWidth', 2);hold on;
plot(th(91, :).*180/pi, db(E_FF(91, :)), 'LineWidth', 2);

C_spillover = 1/(2 .* zeta);

f_pattern_square = abs(Eth).^2 + abs(Eph).^2;

U_feed = C_spillover .* f_pattern_square;

%% 

for i = 1:length(d)
    
    drho = d(i)/100; dphi = pi/ 180;

    [rho, phi] = meshgrid(eps:drho:d(i)/2, eps:dphi:2*pi);
    
    theta_ = 2 * atan(rho/(2 * focal_length));
    
    f_hash = focal_length./d(i);
    
    % Spillover efficiency (eta_s)
    
    theta_0 = 2 * acot(4 * f_hash);
    
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
    
    [Eth_, Eph_] = Feed_FF_Superposition(ModeNumberAper, Gamma, Dm, theta_, phi, F, er, mur, R, Transmission_sum);

   
    
    Erho_ = - Eth_;
    Eph_ = - Eph_;
   
    Int_etp_n_rho = Erho_ .* rho .* drho .* dphi;
    etp_n_rho(i) = sum(sum(Int_etp_n_rho));
    
    Int_etp_n_phi = Eph_ .* rho .* drho .* dphi;
    etp_n_phi(i) = sum(sum(Int_etp_n_phi));
    
    etp_n(i) = abs(etp_n_rho(i)).^2 + abs(etp_n_phi(i)).^2;
    
    E_abs_rpz = abs(Erho_).^2 + abs(Eph_).^2;
    Int_etp_d = E_abs_rpz .* rho .* drho .* dphi;
    etp_d(i) = sum(sum(Int_etp_d));
    
    
    Area(i) = pi .* (d(i)^2)/4;
    
    e_tp(i) =  (1/Area(i)) * etp_n(i)./etp_d(i);
    

end

e_ap = eta_s .* e_tp;

figure;
plot(d, eta_s, 'LineWidth', 2);
hold on;
plot(d, e_tp, 'LineWidth', 2);
hold on;
plot(d, e_ap, 'LineWidth', 2);

grid on;

