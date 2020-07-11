function [Max_Exp_diff] = MinXP_Goal_V2L(R, Len, F, k)
c0 = 3e8;
n = length(R);

l1 = c0/F/4;
L = [l1 ones(1, n-1) * Len/n];

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(eps:dth:pi/2+eps, eps:dph:2*pi+eps);

[ModeNumberAper, Transmission_sum]  = Feed_Gamma_Dm_opt_V2L(R, L, F, k, er, mur);

[~, ~, Eco, Exp, ~, ~] = Feed_FF_Superposition_V2L(ModeNumberAper, th, ph, F, er, mur, R, Transmission_sum);

Max_Exp_diff = -abs(max(max(abs(Eco))) - max(abs(Exp(46, :))));

end