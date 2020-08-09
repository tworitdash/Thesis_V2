
clear;
a = 3.4e-2;
b = 3.4e-2;

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

beta_z10 = -1j .* sqrt(-(k1.^2 - kc1.^2));
Y10 = beta_z10./(omega .* mu1);


M = 5;
L = 100;

for o = 1:length(L)
    disp('Iteration over L');
    disp(o);

for p = 1:length(M)
    
m = 1:1:M(p);

modes = 2 * m + 1;

for q = 1:length(m)
    kci = sqrt((modes(q) * pi./a).^2 + (n_pro .* pi./b).^2);
    beta_z = -1j .* sqrt(-(k1.^2 - kci.^2));
    Ym0(q) = beta_z./(omega.*mu1);
end

for i = 1:length(m)
   for j = 1:length(m)
    if i == j
            
        Ymn_mut(i, j) = Ymn_numerical(k2, a, b, modes(i), modes(j), omega, mu2, L(o)) + Ym0(i);
        
    else
        Ymn_mut(i, j) = Ymn_numerical(k2, a, b, modes(i), modes(j), omega, mu2, L(o));
    end
   end 
end

for l = 1:length(m)
    Yrhs(l) = Ymn_numerical(k2, a, b, 1, modes(l), omega, mu2, L(o));
end

Dm = Ymn_mut\(-Yrhs.');

Dm_(o, :) = Ymn_mut\(-Yrhs.');

Y11(o) = Ymn_numerical(k2, a, b, 1, 1, omega, mu2, L(o));

yap(o) = (Y11(o) + Dm.' * Yrhs.')./Y10;

yap_(o) = (Dm.' * Yrhs.')./Y10;


Gamma(o) = (1 - yap(o))./(1 + yap(o));

end
end

% for p = 1:5
%     hold on;
%     plot(m, db(abs(Dm_(p, :))), 'Linewidth', 2);
%     hold on;
% end
figure;
hold on;
plot(L, db(abs(Gamma)), 'Linewidth', 2);
grid on;
% 
figure;
hold on;
plot(L, angle(Gamma)*180/pi, 'Linewidth', 2);
grid on;

figure;
plot(modes, db(abs(Dm)), 'Linewidth', 2);
grid on;

% figure;
% plot(L, (abs(yap)), 'Linewidth', 2);
% grid on;
% 
% figure;
% plot(L, (abs(yap_)), 'Linewidth', 2);
% grid on;
% 
% figure;
% plot(L, real((yap - yap_).*Y10), 'Linewidth', 2);
% hold on;
% plot(L, imag((yap - yap_).*Y10), 'Linewidth', 2);
% grid on;
