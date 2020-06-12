

F = 5e9;
c0 = 3e8;

omega = 2 * pi * F;

Y0 = 1./(120 * pi);

k0 = omega./c0;

lamb = 2 * pi ./ k0;

x = 0.6:0.01:1;

% R = 2e-2;

R =  x .* lamb ./ 2;

for i = 1:length(R)

    r = R(i);

Str = load('Xmn.mat');
Xmn = Str.Xmn;
xmn = Xmn(1).xmn;


C1 = 2.*(xmn)^2.*(xmn./(k0 .* r))^2./((xmn.^2 - 1).*sqrt(1 - (xmn./(k0 .*r))^2));

C2 = 2./((xmn.^2 - 1).*sqrt(1 - (xmn./(k0 .* r)).^2));

I1 = @(beta) beta .* (-1j.*sqrt(-(1 - beta.^2))) .* (besselj_der(1, k0.*r.*beta)).^2./((xmn./(k0 .* r)).^2 - beta.^2).^2;
I2 = @(beta) (besselj(1, k0 .* r .* beta)).^2./(beta .* (-1j.*sqrt(-(1 - beta.^2))));

del = 0.01.*k0;

wpoints = [eps (1+1j).*eps k0 + 1j.*eps];

func = @(beta) C1 .* I1(beta) + C2 .* I2(beta);

% yin = conj(integral(func, eps, 50.*k0+1j.*eps, 'Waypoints', wpoints));

yin_(i) = integral(func, eps, 1+1j.*eps);

yin_1(i) = integral(func, 1+1j.*eps, 100+1j.*eps);

yin(i) = real(yin_(i)) + 1j.*imag(yin_1(i));

yap(i) = yin(i) .* sqrt(1 - (xmn./(k0.*r))^2);


% Gamma(i) = (1 - yap(i))/(1 + yap(i));
Gamma_(i) = (1 - yin(i))./(1 + yin(i));

end

figure(1);
hold on;
plot(2*R./lamb, real(yin_), 'Linewidth', 2);

hold on;
plot(2*R./lamb, imag(yin_1), 'Linewidth', 2);
grid on;


xlabel('2r/ \lambda', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Y_{in}', 'FontWeight', 'bold', 'FontSize', 16);
title('Input Admittance', 'FontWeight', 'bold', 'FontSize', 16);

legend({'g_{in} K space integrals', 'b_{in} K space integrals', 'Mishustin integrals', 'Mishustin integrals'...
    }, 'location', 'northeast', 'FontWeight', 'bold', 'FontSize', 16);


figure(3);
hold on;
% plot(2*R./lamb, abs(Gamma), 'Linewidth', 2);
% hold on;
plot(2*R./lamb, db(abs(Gamma_)), 'Linewidth', 2);
grid on;

xlabel('2r/ \lambda', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('\Gamma in dB', 'FontWeight', 'bold', 'FontSize', 16);
title('Reflection coefficient', 'FontWeight', 'bold', 'FontSize', 16);

legend({'K Space integrals 19 higher order modes', 'Mishustin integrals only with TE11/TE11'...
    , 'K Space integrals only with TE11/TE11'}, 'location', 'northeast', 'FontWeight', 'bold', 'FontSize', 16);

figure(4);
hold on;
% plot(2*R./lamb, abs(Gamma), 'Linewidth', 2);
% hold on;
plot(2*R./lamb, (angle(Gamma_)).*180/pi, 'Linewidth', 2);
grid on;

xlabel('2r/ \lambda', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Phase of \Gamma in deg ', 'FontWeight', 'bold', 'FontSize', 16);
title('Reflection coefficient', 'FontWeight', 'bold', 'FontSize', 16);

legend({'K Space integrals 19 higher order modes', 'Mishustin integrals only with TE11/TE11'...
    , 'K Space integrals only with TE11/TE11'}, 'location', 'northeast', 'FontWeight', 'bold', 'FontSize', 16);
