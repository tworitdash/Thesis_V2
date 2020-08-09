function [P_der] = legendre_der_(n, m, P_, z)

    P_nminus1 = legendre(n - 1, z);

    P_der = (1./(z^2 - 1)) .* (z .* n .* P_(m, 1) - (m + n) .* P_nminus1(m, 1));
end