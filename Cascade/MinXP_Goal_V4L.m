function [Max_Exp_diff, Eth, Eph, th, ph, Gamma, Dm] = MinXP_Goal_V4L(R_cone, Len, F, k, timesk0, HM)

c0 = 3e8;
lamb = c0./5e9;

% R1 = R_vec(1);
% Rend = R_vec(2);

% R = linspace(R1, Rend, round(Len/(lamb/10)));

% R = linspace(R1, Rend, E);
lamb_opt_freq = c0./F;

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

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(-pi/2:dth:pi/2+eps, eps:dph:2*pi+eps);

[Gamma, Dm, ModeNumberAper, Transmission_sum]  = Feed_Gamma_Dm_opt(R, L, F, k, er, mur, timesk0, HM);
% [ModeNumberAper, Transmission_sum]  = Feed_Gamma_Dm_opt_V2L(R, L, F, k, er, mur);


[Eth, Eph, Eco, Exp, ~, ~] = Feed_FF_Superposition_V3paper(ModeNumberAper, Gamma, Dm, th, ph, F, er, mur, R, Transmission_sum, HM);
% [Eth, Eph, Eco, Exp, ~, ~] = Feed_FF_Superposition_V2L_paper(ModeNumberAper, th, ph, F, er, mur, R, Transmission_sum);



% Max_Exp_diff = -abs(max(max(abs(Eco)))./max(abs(Exp(46, :))));
Max_Exp_diff = db(max(abs(Exp(46, :)))./abs(max(max(abs(Eco)))));

end