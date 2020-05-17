F = 3.3e9;

omega = 2 * pi * F;

mur = 1;
mu0 = 1.25663706e-6;  % Free Space Permeability
mu = mu0 .* mur;

c0 = 3e8;

k0 = omega ./ c0;

a = 3.42e-2;
b = 3.42e-2;

dx = a/200;
dy = b/200;

ep = 0.001*a;

x_i = ep:dx:a;
y_i = ep:dy:b;

[x, y] = meshgrid(x_i, y_i);

m = 1;
n = 1;

C = 1j .* 8 .* a ./ (b .* omega .* mu);

ph = exp(-1j .* k0 .* sqrt(x.^2 + y.^2))./(sqrt(x.^2 + y.^2));

[x_term] = fmn_rec(m, n, a, x, k0);
[y_term] = b - y;

Ymn_i = C .* y_term .* x_term .* ph;

trapz(x_i, trapz(y_i, Ymn_i))

sum(sum(Ymn_i * dx * dy))

figure;
surface(x(1, :)/a, y(:, 1)/b, abs(Ymn_i)); shading flat;
colorbar;
colormap('jet');

% Ymn_mut = sum(sum(Ymn_i));

