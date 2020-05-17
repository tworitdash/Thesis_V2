function [fmn] = fmn_rec(m, n, a, x, k0)
    if m == n
        C = 1./(4 .* pi^2 .* a .* n);
        T1 = 1./pi .* (k0.^2 + (n .* pi./a).^2) .* sin(n .* pi .* abs(x)./a);
        T2 = n./a .* (k0.^2 - (n .* pi./a).^2) .* (a - abs(x)) .* cos(n .* pi .* abs(x)./a);
        fmn = C .* (T1 + T2);
    else
        C = -1j.^(m + n)./(2 * pi^2 * a * (n^2 - m^2));
        T1 = (k0.^2 - (m .* pi ./a).^2) .* n .* sin(m .* pi .* abs(x)./a);
        T2 = (k0.^2 - (n .* pi ./a).^2) .* m .* sin(n .* pi .* abs(x)./a);
        fmn = C .* (T1 - T2);
    end
end