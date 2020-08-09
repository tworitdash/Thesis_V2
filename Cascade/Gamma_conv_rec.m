
    
c0 = 3e8;

a = 3e-2;
b = 1e-2;

mur = 1;

m_pro = 1;
n_pro = 0;

kc = sqrt((m_pro.*pi./a).^2 + (n_pro.*pi./b).^2);

fc = c0 .* kc ./ (2 * pi);

F = fc+1e9;

omega = 2 .* pi .* F;

k0 = omega./c0;

odd = 1:1:30;

m = 2 .* odd + 1;
% m = odd;
n = 0;

[Ym0] = Ymode(m, n, omega, mur, a, b);

N = 1:100:2000;

for o = 1:length(N)

dx = a/N(o);
dy = b/N(o);

x = -a/2-1e-7:dx:a/2-1e-7;
y = -b/2-1e-7:dy:b/2-1e-7;

[x, y] = meshgrid(x, y);

for i = 1:1:length(m)
    for j = 1:1:length(m)
        if i == j
            [Ymn_mut(i, j)] = ymn(m(i), m(j), a, b, x+a/2, y+b/2, k0, mur, dx, dy) + Ym0(i);
        else
            [Ymn_mut(i, j)] = ymn(m(i), m(j), a, b, x+a/2, y+b/2, k0, mur, dx, dy);
        end

    end
end

for l = 1:1:length(m)
    Y_rhs(l) = -ymn(1, m(l), a, b, x+a/2, y+b/2, k0, mur, dx, dy);
end

Dm = Ymn_mut\Y_rhs.';

[Y10, kz10] =  Ymode(1, n, omega, mur, a, b);

[Y11] = ymn(1, 1, a, b, x+a/2, y+b/2, k0, mur, dx, dy);

yap(o) = 1./Y10 .* (Y11 + Dm.' * (-Y_rhs).');

% Gamma = (1 - real(yap))./(1 + real(yap));

Gamma(o) = (1 - (yap(o)))/(1 + (yap(o)));

end

figure;
plot(N, abs(Gamma), 'LineWidth', 2);