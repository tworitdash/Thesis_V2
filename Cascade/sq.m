F = 5e9;
c0 = 3e8;
k0 = 2*pi*F/c0;

krho_re = -2*k0:k0/100:2*k0;
krho_im = -2*k0:k0/100:2*k0;

[krho_re, krho_im] = meshgrid(krho_re , krho_im);

krho = krho_re + 1j .* krho_im;

kz = -1j .* sqrt(-(k0.^2 - krho.^2));
figure;
surf(real(krho)/k0, imag(krho)/k0, real(kz));shading flat; colormap('jet'); colorbar;

xlabel('Reak(k_{\rho})/k_0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Imag(k_{\rho})/k_0', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Real(k_z)', 'FontSize', 12, 'FontWeight', 'bold');


figure;
surf(real(krho)/k0, imag(krho)/k0, imag(kz));shading flat; colormap('jet'); colorbar;

xlabel('Reak(k_{\rho})/k_0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Imag(k_{\rho})/k_0', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Imag(k_z)', 'FontSize', 12, 'FontWeight', 'bold');