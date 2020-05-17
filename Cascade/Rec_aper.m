clear;
c0 = 3e8;


er = 2.53; % Rexolite
er2 = 1; % Freespace

v = c0./sqrt(er); %Rexolite
v2 = c0./sqrt(er2); % Freespace

a = 3.42e-2;
b = 3.42e-2;

mur = 1;

m_pro = 1;
n_pro = 0;

kc = sqrt((m_pro.*pi./a).^2 + (n_pro.*pi./b).^2);

fc = v .* kc ./ (2 * pi);

F = fc+1e9;
% F = 3.3e9;

omega = 2 .* pi .* F;

kr = omega./v; % Rexolite
k0 = omega./c0;

odd = 1:1:15;

m = 2 .* odd + 1;
% m = odd;
n = 0;

[Ym0, ~] = Ymode(m, n, omega, mur, er, a, b);

dx = a/2000;
dy = b/2000;

ep = 1e-7;

x = -a/2-ep:dx:a/2-eps;
y = -a/2-ep:dy:b/2-eps;

[x, y] = meshgrid(x, y);

for i = 1:1:length(m)
    for j = 1:1:length(m)
        if i == j
            [Ymn_mut(i, j)] = ymn(m(i), m(j), a, b, x+a/2, y+b/2, k0, mur, er2, dx, dy) + Ym0(i);
%               [Ymn_mut(i, j)] = ymn_k(m(i), m(j), a, b, k0, mur, er) + Ym0(i);
        else
            [Ymn_mut(i, j)] = ymn(m(i), m(j), a, b, x+a/2, y+b/2, k0, mur, er2, dx, dy);
%             [Ymn_mut(i, j)] = ymn_k(m(i), m(j), a, b, k0, mur, er);
        end

    end
end

for l = 1:1:length(m)
    Y_rhs(l) = -ymn(1, m(l), a, b, x+a/2, y+b/2, k0, mur, er2, dx, dy);
%       Y_rhs(l) = -ymn_k(1, m(l), a, b, k0, mur, er);
end

Dm = Ymn_mut\Y_rhs.';

[Y10, kz10] =  Ymode(1, n, omega, mur, er, a, b);

[Y11] = ymn(1, 1, a, b, x+a/2, y+b/2, k0, mur, er2, dx, dy);
% Y11 = ymn_k(1, 1, a, b, k0, mur, er);

yap = 1./Y10 .* (Y11 + Dm.' * (-Y_rhs).');

% Gamma = (1 - real(yap))./(1 + real(yap));

Gamma = (1 - (yap))/(1 + (yap));

ZTE10 = 1./Y10;

E0y = (-1j .* kz10)./kc.^2 .* (m_pro * pi./a) .* ZTE10;

% EyTE = E0y .* (1 + Gamma) .* cos(m_pro .* pi .* x./a) + sum(E0y .* Dm .* (1 + Gamma) * cos(m .* pi .* x./a));

EyTE_fundamental = E0y .* cos(m_pro .* pi .* x./a);
EyTE = EyTE_fundamental .* (1 + (Gamma));

for o = 1:length(m)
    m_o = m(o);
    n_o = 0;
    [Yo0, kzo0] =  Ymode(m_o, n_o, omega, mur, a, b);
    Zo0 = 1./Yo0;
    kco = sqrt((m_o.*pi./a).^2 + (n_o.*pi./b).^2);
    Eoy = (-1j .* kzo0)./kco.^2 .* (m_o * pi./a) .* Zo0;
    
    EyTE = EyTE + Eoy .* Dm(o) .* (1 + (Gamma)) * cos(m_o .* pi .* x./a);
end

% figure;
% 
% surface(x, y, db(abs(EyTE_fundamental)./max(max(abs(EyTE_fundamental))))); shading flat;
% colormap('jet');
% 
% figure;
% surface(x, y, db(abs(EyTE)./max(max(abs(EyTE))))); shading flat;
% colormap('jet');

figure;
surface(x, y, db(abs(EyTE))); shading flat;
colormap('jet');
colorbar;

% figure;
% 
% surface(x, y, db(abs(EyTE_fundamental))); shading flat;
% colormap('jet');
% colorbar;

% figure;
% 
% surface(x, y, (abs(EyTE_fundamental) - abs(EyTE))); shading flat;
% colormap('jet');
% colorbar;

nomfic = '/Users/tworitdash/course/feko/NF_recv2_FS_MF_lam10.efe';

[EF_xFS, EF_yFS, EF_zFS, x_f, y_f] = E_FEKO_rec(nomfic);
% 
% nomfic = '/Users/tworitdash/course/feko/NF_recv2_abs.efe';
% 
% [EF_xabs, EF_yabs, EF_zabs, x_f, y_f] = E_FEKO_rec(nomfic);

figure;
surface(x_f, y_f, db(abs(EF_yFS.'))); shading flat;
colormap('jet');
colorbar;

% figure;
% 
% surface(x_f, y_f, db(abs(EF_yabs.'))); shading flat;
% colormap('jet');
% colorbar;

% figure;
% 
% surface(x_f, y_f, (abs(EyTE_fundamental(1:1000, 1:1000)) - abs(EF_yabs).')); shading flat;
% colormap('jet');
% colorbar;

figure;

surface(x_f, y_f, db(abs(EyTE(1:1000, 1:1000)) - abs(EF_yFS).')); shading flat;
colormap('jet');
colorbar;

%% Radiation pattern
dth = pi/180;
dphi = pi/180;

z = 5e-2;

ExTE = zeros(size(x));
EzTE = zeros(size(x));

[theta, phi] = meshgrid(-pi/2-eps:dth:pi/2-eps, eps:dphi:2*pi);

[Eth, Eph] = FF_M_cart(E0y, theta, phi, F, m, Dm, Gamma, a, b);

E_abs = sqrt(abs(Eth).^2 + abs(Eph).^2);

[Eth_absorb, Eph_absorb] = FF_M_cart(E0y, theta, phi, F, 0, 0, 0, a, b);
% [Eth_absorb, Eph_absorb] = FF_M_cart_numerical(ExTE, EyTE_fundamental, EzTE, x, y, z, theta, phi, F, dx, dy);

E_abs_absorb = sqrt(abs(Eth_absorb).^2 + abs(Eph_absorb).^2);

[EFth_FS, EFph_FS] = FF_M_cart_numerical(EF_xFS.', EF_yFS.', EF_zFS.', x_f, y_f, z, theta, phi, F, dx, dy);
EF_FSabs = sqrt(abs(EFth_FS).^2 + abs(EFph_FS).^2);

[EFth_abs, EFph_abs] = FF_M_cart_numerical(EF_xabs.', EF_yabs.', EF_zabs.', x_f, y_f, z, theta, phi, F, dx, dy);
EF_absabs = sqrt(abs(EFth_abs).^2 + abs(EFph_abs).^2);

% patternCustom(E_abs,theta,phi);

%FF Plots


figure(100);

hold on;  
% plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))/max(abs(E_abs(1, :)))), 'LineWidth', 2);
plot(theta(1, :)*(180/pi), db(abs(E_abs_absorb(1, :))), 'LineWidth', 2);
hold on;
plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))), '-.', 'LineWidth', 2);
hold on;
plot(theta(1, :)*(180/pi), db(abs(EF_absabs(1, :))), '*', 'LineWidth', 0.5);
hold on;
plot(theta(1, :)*(180/pi), db(abs(EF_FSabs(1, :))), '--', 'LineWidth', 1);
hold on;

% figure(42)
hold on;
% plot(theta(91, :)*(180/pi), db(abs(E_abs(91, :))/max(abs(E_abs(91, :)))), 'LineWidth', 2);
plot(theta(91, :)*(180/pi), db(abs(E_abs_absorb(91, :))), 'LineWidth', 2);
hold on;
plot(theta(91, :)*(180/pi), db(abs(E_abs(91, :))), '-.', 'LineWidth', 2);
hold on;
plot(theta(91, :)*(180/pi), db(abs(EF_absabs(91, :))), '*', 'LineWidth', 0.5);
hold on;
plot(theta(91, :)*(180/pi), db(abs(EF_FSabs(91, :))), '--', 'LineWidth', 1);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');


grid on;

ylim([-50 20]);
xlim([-90 90]);


legend;

