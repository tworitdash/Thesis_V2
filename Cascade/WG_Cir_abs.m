%% Far fields of the ctlindrical waveguide analytically for TE11 mode
clear;
c0 = 3e8;

er = 1; mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;



dth = pi/180;
dphi = pi/180;

F = 5e9;
beta = 2 * pi * F / c0;

[theta, phi] = meshgrid(-pi/2-eps:dth:pi/2-eps, eps:dphi:2*pi);

R = 2e-2;

Str = load('Xmn.mat');

str = Str.Xmn;

xmn = str(1).xmn;
m = str(1).m;
mode = str(1).mode;

beta_rhop = xmn./R;

ZTE = 2 * pi * F * mu ./ beta;

% if m == 0
%   deltam = 1;
% else
%   deltam = 0;
% end
% if mode == "TE"
%    Nup = (pi*(1+deltam)/2 .* (xmn.^2 - m.^2) .* (besselj(m, xmn)).^2).^(-1);
% elseif modep  == "TM"
%    Nup = (pi*(1+deltam)/2 .* (xmn_p).^2 .* (besselj_der(m, xmn)).^2).^(-1);
% end
Nup = 1j .* beta ./ beta_rhop.^2 * ZTE;

k0 = (2 * pi * F)/c0;

kx = k0 .* sin(theta) .* cos(phi);
ky = k0 .* sin(theta) .* sin(phi);

kz = (-1j)*sqrt(-(k0.^2 - kx.^2 - ky.^2));

Nu = sqrt((kx.^2 + ky.^2));

Psi = acos(kx./Nu);

for i = 1:size(Nu, 1)
    
    for k = 1:size(Nu, 2)

        I0(i, k) = Lommel(0, R, beta_rhop, Nu(i, k), m - 1, m - 1);
        I2(i, k) = Lommel(0, R, beta_rhop, Nu(i, k), m + 1, m + 1);

    end
end

E_ft_x = (Nup) .* pi .* beta_rhop .* sin(2 * Psi) .* I2;
E_ft_y = (Nup) .* pi .* beta_rhop .* (I0 + cos(2 .* Psi) .* I2);


c2 = 1j .* k0 / (4 * pi);

Eth = c2 .* (E_ft_x .* cos(phi) + E_ft_y .* sin(phi));
Eph = c2 .* cos(theta) .* (E_ft_y .* cos(phi) - E_ft_x .* sin(phi));


E_abs = sqrt(abs(Eth).^2 + abs(Eph).^2);

figure(96);

hold on;
plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))), '-.', 'LineWidth', 2);

plot(theta(91, :)*(180/pi), db(abs(E_abs(91, :))), '-.', 'LineWidth', 2);



xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');


grid on;

ylim([-50 10]);
xlim([-90 90]);


legend;

