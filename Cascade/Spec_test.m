
k0 = 60;

dkx = k0/100;

kx = -5*k0:dkx:5*k0;

b = 2e-2;
a = 5e-2;

for i = 1:length(kx)
    
kxi = kx(i);

del = 0.01 .* k0;
    
% dy = k0/100;

% ky = - 5 * k0:dky:5 * k0;

Y = @(ky) (k0.^2 - kxi.^2)./(-1j .* sqrt(-(k0.^2 - kxi.^2 - ky.^2))) .* (besselj(0, ky.*b./2));

% Y = @(ky) (k0^2 - kxi^2)./(-1j .* sqrt(-(k0.^2 - kxi.^2 - ky.^2))) .* (sinc(ky .* b./2/pi)).^2;

y(i) = integral(Y, -50.*k0-1j*del,50.*k0+1j*del, 'Waypoints', [(-1-1j).*del, (1+1j).*del]);

end

% Int = @(kxy) Dkx .* cos(kxy.*a/2) ./ (((m * pi).^2 - (kxy .* a).^2) .* ((n * pi).^2 - (kxy .* a).^2));
% 
% Integr = integral(Int, -50.*k0-1j*del,50.*k0+1j*del, 'Waypoints', [(-1-1j).*del, (1+1j).*del]);

y_analyt = pi .* (k0.^2 - kx.^2) .* besselj(0, b./4 .* -1j .* sqrt(-(k0.^2 - kx.^2))) .* besselh(0, 2,  b./4 .* -1j .* sqrt(-(k0.^2 - kx.^2))); 
% figure;
% plot(ky./k0, abs(Y), 'Linewidth', 2);
% grid on;
figure;
hold on;
plot(kx./k0, abs(y), 'Linewidth', 2);
hold on;
grid on;
plot(kx./k0, abs(y_analyt), '-.', 'Linewidth', 2);