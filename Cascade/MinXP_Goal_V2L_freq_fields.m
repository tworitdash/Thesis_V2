function [Eth_, Eph_, Eco_, Exp_, CO_, XP_, E, th, ph, Max_Exp_diff] = MinXP_Goal_V2L_freq_fields(R_cone, Len, F, k)

c0 = 3e8;
lamb = c0./5e9;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
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

% c0 = 3e8;
% n = length(R);
% 
% l1 = c0/F/4;
% L = [l1 ones(1, n-1) * Len/n];

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG section

epsilon = er .* er0;
mu = mur .* mu0;

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(-pi/2+eps:dth:pi/2+eps, eps:dph:2*pi+eps);

[ModeNumberAper, Transmission_sum]  = Feed_Gamma_Dm_opt_V2L_freq(R, L, F, k, er, mur);

[rho, phi] = meshgrid(eps:R(end)/100:R(end), eps:pi/180:2*pi);

parfor i = 1:length(F)
    [Er_rho, Er_phi, ~] = E_n(1:1:ModeNumberAper, rho, phi, F(i), R(end), L(end), epsilon(end), mu(end));
    T = Transmission_sum(i).str_req;
    
    Eth = zeros(size(th));
    Eph = zeros(size(ph));
    for k = 1:ModeNumberAper
        [Eth1, Eph1] = FF_M_V2(squeeze(Er_rho(k, :, :)), squeeze(Er_phi(k, :, :)), rho, phi, th, ph, F(i), R(end)/100, pi/180);
        Eth = Eth + T(k) .* Eth1;
        Eph = Eph + T(k) .* Eph1; 
    end
    
    Exp = cos(ph) .* Eth - sin(ph) .* Eph;
    Eco = sin(ph) .* Eth + cos(ph) .* Eph;
    
%     [Eth, Eph, Eco, Exp, CO, XP] = Feed_FF_Superposition_V2L(ModeNumberAper(i, end), th, ph, F(i), er, mur, R, Transmission_sum(i).str_req);

% Max_Exp_diff = -abs(max(max(abs(Eco)))./max(abs(Exp(46, :))));
%     Max_Exp_diff(i) = db(max(abs(Exp(46, :)))./abs(max(max(abs(Eco)))));
    Eth_(i, :, :) = Eth;
    Eph_(i, :, :) = Eph;
    Eco_(i, :, :) = Eco;
    Exp_(i, :, :) = Exp;
%     CO_(i, :, :) = CO;
%     XP_(i, :, :) = XP;
    CO_(i, :, :) = zeros(size(th));
    XP_(i, :, :) = zeros(size(ph));
    E(i, :, :) = sqrt(abs(Eth).^2 + abs(Eph).^2);
    Max_Exp_diff(i) = db(max(abs(Exp(46, :)))./abs(max(max(abs(Eco)))));
end

end