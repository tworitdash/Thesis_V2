c0 = 3e8;

F = 21e9; % Frequency at which far field is requested
lamb = c0/F;



rt = 2 * lamb; % radius

Length = 4 * lamb;

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

% for i = 1:n
%     f =  fc(R(i), er(i), mur(i));
%     N_i  =  find(f < F);
%     N(i) = length(N_i);
% end
N = 1;

% L = ones(1, n) .* Length/n; % length of each waveguide section
% L(1) = 5 * L(1);
% L(end) = 7.5e-3;

z = L;

[Ep_rho, Ep_phi] = E_n(1:1:N(end), rho, ph, F, R(end), z, epsilon(end), mu(end), drho, dphi);


E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
% E_aperture_z = zeros(size(rho));

for k = 1:N(end)
    
%      if (k == 2) || (k == 10)
%         j = 1;
%     else
%         j = 2;
%     end
%     if (k == 4)
%         l = 2;
%     else
%         l = 1;
%     end
    E_aperture_rho = E_aperture_rho + squeeze(Ep_rho(k, :, :))  .* exp(1j .* pi/2);% .* A(l, k);
    E_aperture_phi = E_aperture_phi + squeeze(Ep_phi(k, :, :))  .* exp(1j .* pi/2);% .* A(j, k);
%     E_aperture_z = E_aperture_z + squeeze(Ep_z(k, :, :)) .* (Gamma_sum(k));
end



%% FF
[theta, phi] = meshgrid(-pi/2-eps:pi/180:pi/2-eps, eps:pi/180:2*pi+eps);
    


[Eth, Eph] = FF_M_V2(E_aperture_rho, E_aperture_phi, rho, ph, theta, phi, F, drho, dphi);


E_abs = sqrt(abs(Eth).^2 + abs(Eph).^2);

figure(41);

hold on;
plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))/max(abs(E_abs(1, :)))), 'LineWidth', 2);
hold on;
plot(theta(91, :)*(180/pi), db(abs(E_abs(91, :))/max(abs(E_abs(91, :)))), 'LineWidth', 2);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');

grid on;

ylim([-50 0]);

% figure(9);
% 
% hold on;
%     
% plot(theta(1, :)*(180/pi), db(abs(E_abs(1, :))/max(abs(E_abs(1, :)))), 'LineWidth', 2);
% hold on;
% plot(theta(90, :)*(180/pi), db(abs(E_abs(90, :))/max(abs(E_abs(90, :)))), 'LineWidth', 2);
% 
% 
% xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('E_{abs} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Far electric field', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'\phi = 0', '\phi = 90'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
% grid on;
% 
% ylim([-50 0]);
