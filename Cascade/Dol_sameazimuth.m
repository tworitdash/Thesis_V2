function [r, dol] = Dol_sameazimuth(f, er, mur) % D over lambda
 % er = 1;
% mur = 1;
c0 = 3e8;
lamb = c0/f;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium


mu = mu0 * mur;

Str = load('Xmn_azimuthal_inc_TE.mat');
Xmn = Str.xmn_TE;

for i = 1:size(Xmn(1:318), 2)
    X(i) = Xmn(i).xmn;
    M(1, i) = Xmn(i).mode;
    M(2, i) = Xmn(i).m;
    M(3, i) = Xmn(i).n;
end


r = X ./ (2 * pi * f * sqrt(mu .* epsilon));

dol = (2 .* r)/lamb;
end