function [Ymn_mut] = ymn_k(m, n, a, b, k0, mur, er)

c0 = 3e8;
omega = c0 .* k0;

mu0 = 1.25663706e-6;  % Free Space Permeability
mu = mu0 .* mur;

dkx = k0/100;
dky = k0/100;

kx = cat(2, -20*k0:dkx:-k0, k0:dkx:20*k0);
ky = cat(2, -20*k0:dky:-k0, k0:dky:20*k0);
% ky = -10*k0:dky:10*k0;

[kx, ky] = meshgrid(kx, ky);

kz = -1j .* sqrt(-(k0.^2 - kx.^2 - ky.^2));

C = -1j^(m + n) * m .* n;
N = 2 .* a .* b .* (k0.^2 - kx.^2) .* (sin(ky * b ./ 2)).^2 .* (cos(kx .* a ./2)).^2 .* dkx .* dky;
D = omega .* mu .* kz .* (ky .* b./2).^2 .* ((m .* pi).^2 - (kx .* a).^2) .* ((n .* pi).^2 - (kx .* a).^2);



Ymn_i = C .* N./D;

Ymn_mut = sum(sum(Ymn_i));

% figure;
% surf(kx(1, :)./k0, ky(:, 1)./k0, abs(Ymn_i)); shading flat;
% colorbar;
% colormap('jet');

end
