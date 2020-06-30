function [Max_Exp] = MinXP_Goal(R, L, F, k, timesk0, HM)

n = length(R);


er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(eps:dth:pi/2+eps, eps:dph:2*pi+eps);

[Gamma, Dm, ModeNumberAper, Transmission_sum]  = Feed_Gamma_Dm_opt(R, L, F, k, er, mur, timesk0, HM);

[~, ~, ~, Exp, ~, ~] = Feed_FF_Superposition(ModeNumberAper, Gamma, Dm, th, ph, F, er, mur, R, Transmission_sum);

Max_Exp = max(Exp(46, :));

end