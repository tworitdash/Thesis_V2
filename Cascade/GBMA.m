
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


Eax = cos(ph) .* Ep_rho - sin(ph) .* Ep_phi;
Eay = sin(ph) .* Ep_rho + cos(ph) .* E_phi;

%% Free space GBMA

Psi_0 = Psi(0, R, L, F, rho, phi, z_obs); 
Psi_2 = Psi(2, R, L, F, rho, phi, z_obs);

