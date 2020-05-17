
k0 = 60;

dkx = k0/100;

kx = -5*k0:dkx:5*k0;

b = 2e-2;
a = 5e-2;

del = 0.01 .* k0;
    
% dy = k0/100;

% ky = - 5 * k0:dky:5 * k0;

% Y = @(ky) (k0.^2 - kxi.^2)./(-1j .* sqrt(-(k0.^2 - kxi.^2 - ky.^2))) .* (besselj(0, ky.*b./2));
m = 1;
n = 3;

Y = @(kxy) Dkx(kxy, k0, b) .* cos(kxy.*a/2) ./ (((m * pi).^2 - (kxy .* a).^2) .* ((n * pi).^2 - (kxy .* a).^2));

Ymn = integral(Y, -50.*k0-1j*del,50.*k0+1j*del, 'Waypoints', [(-1-1j).*del, (1+1j).*del]);



