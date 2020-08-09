function [Ymn] = Ymn_numerical(k0, a, b, m, n, omega, mu, L)

del = 0.01 .* k0;

disp(m*pi/a);
disp(n*pi/a);

wpoints = [-eps-1j.*eps eps+1j.*eps k0+eps+1j.*eps k0+eps];


Y = @(kxy) Dkx(kxy, k0, b) .* (k0.^2 - kxy.^2) .* (cos(kxy.*a/2)).^2 ./ (((m * pi).^2 - (kxy .* a).^2) .* ((n * pi).^2 - (kxy .* a).^2));

C = -1j^(m + n) .* m .* n .* 2 .* a .* b ./ (omega .* mu);

Ymn = C .* integral(Y, -L.*k0-1j*del,L.*k0+1j*del, 'Waypoints', [(-1-1j).*del, (1+1j).*del]);

% Ymn = C .* integral(Y, -eps-1j.*L.*k0,eps-1j.*L.*k0, 'Waypoints', wpoints);


end


