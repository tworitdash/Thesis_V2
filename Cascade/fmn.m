function  TermX = fmn(x, m, n, a, k)
    if m == n
        C = 1./(4 .* pi.^2 .* a .* n);
        K_x = (1./pi) .* (k.^2 + (n .* pi./a).^2) .* sin(n .* pi .* abs(x)./a);
        K_y = (n./a) .* (k.^2 - (n .* pi./a).^2) .* (a - abs(x)) .* cos(n .* pi .* abs(x)./a);
        TermX = C .* (K_x + K_y);
    else
        C1 = -1j.^(m + n)./(2 .* pi.^2 .* a .* (n^2 - m^2));
        K_x1 = (k.^2 - (m .* pi./a).^2) .* n .* sin(m .* pi .* abs(x)./a);
        K_y1 = (k.^2 - (n .* pi./a).^2) .* m .* sin(n .* pi .* abs(x)./a);
       TermX = C1 .* (K_x1 - K_y1);
    end
end