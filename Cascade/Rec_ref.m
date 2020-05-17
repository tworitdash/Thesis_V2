
clear;
a = 8.3e-2;
b = 1e-2;

c0 = 3e8;

F = 3.3e9;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

er1 = 2.53;
mur1 = 1;

epsilon1 = er1 .* er0;
mu1 = mur1 .* mu0;

er2 = 1;
mur2 = 1;

epsilon2 = er2 .* er0;
mu2 = mur2 .* mu0;

v1 = c0./sqrt(er1);
v2 = c0./sqrt(er2);


omega = 2 * pi * F;

k1 = omega/v1;
k2 = omega/v2;

m_pro = 1;
n_pro = 0;

kc1 = sqrt((m_pro * pi./a).^2 + (n_pro .* pi./b).^2);

fc1 = v1 .* kc1 ./ (2 .* pi);

% ep = 1e-5;

dx = a/1000;
dy = b/1000;

x = eps:dx:a+eps;
y = eps:dy:b+eps;

[x, y] = meshgrid(x, y);

ep = 1e-6;

x__ = ep:dx:a;
y__ = ep:dy:b;
            
[x__, y__] = meshgrid(x__, y__);

beta_z10 = -1j .* sqrt(-(k1.^2 - kc1.^2));
Y10 = beta_z10./(omega .* mu1);

M = 1:1:15;

for p = 1:length(M)
    
m = 1:1:M(p);

modes = 2 * m + 1;

for i = 1:length(m)
    kci = sqrt((modes(i) * pi./a).^2 + (n_pro .* pi./b).^2);
    beta_z = -1j .* sqrt(-(k1.^2 - kci.^2));
    Ym0(i) = beta_z./(omega.*mu1);
end

for i = 1:length(m)
   for j = 1:length(m)
    if i == j
            
            
            Y1(i, j) = Y_mut2(modes(i), modes(j), ep, ep, k2, omega, mu2, epsilon2);
            Y2(i, j) = Y_mut(modes(i), modes(j), x__, y__, a, b, dx, dy, k2, omega, mu2, epsilon2);
           
            Ymn_mut(i, j) = Y1(i, j) + Y2(i, j) + Ym0(i);
%             Ymn_mut(i, j) = Y_mut(modes(i), modes(j), x+a./2, y+b./2, a, b, dx, dy, k2, omega, mu2, epsilon2) + Ym0(i);
        
        
        %           Ymn_mut(i, j) = Y_mut2(modes(i), modes(j), a, b, k2, omega, mu2, epsilon2) + Ym0(i);
    else
           Ymn_mut(i, j) = Y_mut(modes(i), modes(j), x, y, a, b, dx, dy, k2, omega, mu2, epsilon2);
%         Ymn_mut(i, j) = Y_mut2(modes(i), modes(j), a, b, k2, omega, mu2, epsilon2);
    end
   end 
end

for l = 1:length(m)
    
        Yrhs(l) = Y_mut(1, modes(l), x, y, a, b, dx, dy, k2, omega, mu2, epsilon2);
%       Yrhs(l) = Y_mut2(1, modes(l), a, b, k2, omega, mu2, epsilon2);
end

Dm = Ymn_mut\(-Yrhs.');

% Y11 = Y_mut(1, 1, x+a./2, y+b./2, a, b, dx, dy, k2, omega, mu2, epsilon2);

Y111 = Y_mut2(1, 1, ep, ep, k2, omega, mu2, epsilon2);
Y112 = Y_mut(1, 1, x__, y__, a, b, dx, dy, k2, omega, mu2, epsilon2);

Y11 = Y111 + Y112;

yap = (Y11 + Dm.' * Yrhs.')./Y10;


Gamma(p) = (1 - yap)./(1 + yap);

end
figure(2);
hold on;
plot(M, abs(Gamma), 'Linewidth', 2);
grid on;
% 
figure(3);
hold on;
plot(M, angle(Gamma)*180/pi, 'Linewidth', 2);
grid on;
% figure;



%% 


ZTE_10 = 1./Y10;

% %%
% 
% ep = 1e-2;
% 
% dx = a/1000;
% dy = b/1000;
% 
% x = -a/2+ep:dx:a/2+ep;
% y = -b/2+ep:dy:b/2+ep;
% 
% [x, y] = meshgrid(x, y);
% 
% modes = 2 * (1:30) + 1;
% for i = 1:10
% %     for j = 1:15
% %       Y1(i) = Y_mut(5, modes(i), x+a./2, y+b./2, a, b, dx, dy, k2, omega, mu2, epsilon2);% + Ym0(i);
%       Y2(i) = Y_mut2(5, modes(i), a, b, k2, omega, mu2, epsilon2);
% %     end
% end
% 
% % figure;
% % plot(modes, real(Y1)); 
% % hold on;
% % plot(modes, real(Y2));
% % 
% % figure;
% % plot(modes, abs(real(Y1) - real(Y2))./real(Y2) .* 100);
% % 
% % figure;
% % plot(modes, imag(Y1)); 
% % hold on;
% % plot(modes, imag(Y2));
% % 
% % figure;
% % plot(modes, abs(imag(Y1) - imag(Y2))./imag(Y2) .* 100);
% 
% 
% % figure;
% % surface(modes, modes.', real(Y1)); shading flat;
% % colorbar;
% % figure;
% % surface(modes, modes.', imag(Y1)); shading flat;
% % colorbar;
% % 
% % figure;
% % surface(modes, modes.', real(Y2)); shading flat;
% % colorbar;
% % figure;
% % surface(modes, modes.', imag(Y2)); shading flat;
% % colorbar;
