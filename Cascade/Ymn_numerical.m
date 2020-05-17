function [Ymn] = Ymn_numerical(k0, a, b, m, n, omega, mu, L)

del = 0.01 .* k0;


Y = @(kxy) Dkx(kxy, k0, b) .* cos(kxy.*a/2) ./ (((m * pi).^2 - (kxy .* a).^2) .* ((n * pi).^2 - (kxy .* a).^2));

C = -1j^(m + n) .* m .* n .* 2 .* a .* b ./ (omega .* mu);

Ymn = C .* integral(Y, -L.*k0-1j*del,L.*k0+1j*del, 'Waypoints', [(-1-1j).*del, (1+1j).*del]);

end


