c0 = 3e8;

F = 5e9; % Frequency at which far field is requested

rt = 2e-2; % redius

Length = 7e-2;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

R = rt;

drho = R(end)/100;
dphi = pi/180;

[rho, ph] = meshgrid(eps:drho:R(end), eps:dphi:2*pi+eps);

n = length(R);

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

epsilon = er .* er0;
mu = mur .* mu0;

L = Length;

for i = 1:n
    f =  fc(R(i), er(i), mur(i));
    N_i  =  find(f < F);
    N(i) = length(N_i);
end


% L = ones(1, n) .* Length/n; % length of each waveguide section
% L(1) = 5 * L(1);
% L(end) = 7.5e-3;

z = L;

[Ep_rho, Ep_phi, Ep_z] = E_r(1:1:N(end), rho, ph, F, R(end), z, epsilon(end), mu(end));


E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
% E_aperture_z = zeros(size(rho));

for k = 1:N(end)
    E_aperture_rho = E_aperture_rho + squeeze(Ep_rho(k, :, :))  .* exp(1j .* pi/2);
    E_aperture_phi = E_aperture_phi + squeeze(Ep_phi(k, :, :))  .* exp(1j .* pi/2);
%     E_aperture_z = E_aperture_z + squeeze(Ep_z(k, :, :)) .* (Gamma_sum(k));
end



%% FF
[theta, phi] = meshgrid(eps:pi/180:2*pi+eps, eps:pi/180:2*pi+eps);
    


[Ex, Ey, Ez, Hx, Hy, Hz] = FF_M(E_aperture_rho, E_aperture_phi, rho, ph, theta, phi, F, drho, dphi);

H_abs = sqrt(Hx.^2 + Hy.^2 + Hz.^2);

E_abs = sqrt(Ex.^2 + Ey.^2 + Ez.^2);

figure;

    
plot(theta(1, :)*(180/pi), db(abs(H_abs(1, :))/max(abs(H_abs(1, :)))), 'LineWidth', 2);
hold on;
plot(theta(90, :)*(180/pi), db(abs(H_abs(90, :))/max(abs(H_abs(90, :)))), 'LineWidth', 2);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('H_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far magnetic field', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;

ylim([-50 0]);

figure(9);

hold on;
    
plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))/max(abs(E_abs(1, :)))), 'LineWidth', 2);
hold on;
plot(theta(90, :)*(180/pi), db(abs(E_abs(90, :))/max(abs(E_abs(90, :)))), 'LineWidth', 2);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;

ylim([-50 0]);
