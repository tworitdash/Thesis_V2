m = 5; n = 7;

k0 = 33.51;

kx_real = linspace(-2*m*pi/a, 2*m*pi/a, 1000);
kx_imag = linspace(-2*m*pi/a, 2*m*pi/a, 1000);

[kx_real, kx_imag] = meshgrid(kx_real, kx_imag);

kx = kx_real + 1j .* kx_imag;


a = 2e-2;

X = (cos(kx.*a/2)).^2 ./ (((m * pi).^2 - (kx .* a).^2) .* ((n * pi).^2 - (kx .* a).^2));

Y = -1j .* sqrt(-(k0.^2 - kx.^2));

figure;
surface(kx_real(1, :)/k0, kx_imag(:, 1)/k0, real(X)); shading flat;
colormap('jet');
figure;
surface(kx_real(1, :)/k0, kx_imag(:, 1)/k0, imag(X)); shading flat;
colormap('jet');

figure;
surface(kx_real(1, :)/k0, kx_imag(:, 1)/k0, real(Y)); shading flat;
colormap('jet');
figure;
surface(kx_real(1, :)/k0, kx_imag(:, 1)/k0, imag(Y)); shading flat;
colormap('jet');