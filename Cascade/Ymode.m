function [y, kz] = Ymode(M, N, omega, mur, er, a, b)
c0 = 3e8;

v = c0./sqrt(er);

k0 = omega./v;

mu0 = 1.25663706e-6;  % Free Space Permeability
mu = mu0 .* mur;
% ymode = zeros(length(M), length(N));

for i = 1:length(M)
    for j = 1:length(N)
        m = M(i);
        n = N(j);
        kc = sqrt((m.*pi./a).^2 + (n.*pi./b).^2);
        kz(i, j) = -1j .* sqrt(-(k0.^2 - kc.^2));

        y(i, j) = kz(i, j)./(omega .* mu);
    end
end

end