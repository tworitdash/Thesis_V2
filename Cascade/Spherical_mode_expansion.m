clear;
c0 = 3e8;

F = 21e9; % Frequency at which far field is requested
lamb = c0/F;
beta = 2 * pi / lamb;



rt = 2 * lamb; % radius
% rt = 4e-2;

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

%%  --------------------------------------------------------------------------------------------------------

z = 0;

%[Ep_rho, Ep_phi, Ez, Hp_rho, Hp_phi, Hz] = E_r(1:1:N(end), rho, ph, F, R(end), z, epsilon(end), mu(end), drho, dphi);
[Ep_rho, Ep_phi, Ez, Hp_rho, Hp_phi, Hz] = E_r(1:1:N(end), rho, ph, F, R(end), z, epsilon(end), mu(end));


Q_wg = Pow(N, Ep_rho, Ep_phi, Hp_rho, Hp_phi, rho, drho, dphi);


%% Free space spherical modes

r = sqrt(rho.^2 + z^2);
theta = pi/2;

M = 10;
Z0 = 120 * pi;

for i = 1:1:M
    
        S = load('Sper.mat');
        Sper = S.Sper;
        
        mode = Sper(i).mode;
        pol = Sper(i).pol;
        n = Sper(i).n;
        m = Sper(i).m;
    
        H = besselh(n, 2, beta .* r);
        P = legendre(n, cos(theta));
        P_der = legendre_der_(n, m, P, cos(theta));
        
        if mode == "TE"
            if pol == "even"
                Er = zeros(size(rho)) + eps;
                Eth = H .* m .* P(m, 1) ./ sin(theta) .* sin(m .* ph);
                Eph = sin(theta) .* H .* m .* P_der .* cos(m .* ph);
            else
                Er = zeros(size(rho)) + eps;
                Eth = H .* m .* P(m, 1) ./ sin(theta) .* cos(m .* ph);
                Eph = sin(theta) .* H .* m .* P_der .* sin(m .* ph);
            end
        end
        
        
        if mode == "TM"
            if pol == "even"
                Er = n .* (n + 1) .* H ./ (beta .* r) .* P(m, 1) .* sin(m .* ph);
                Eth = - sin(theta) .* 1./ (beta .* r) .* (H + beta .* H_der(n, beta .* r)) .* P_der .* sin(m .* ph);
                Eph = - 1./ (beta .* r) .* (H + beta .* H_der(n, beta .* r)) .* m .* P(m, 1)./sin(theta) .* cos(m .* ph);
            else
                Er = n .* (n + 1) .* H ./ (beta .* r) .* P(m, 1) .* cos(m .* ph);
                Eth = - sin(theta) .* 1./ (beta .* r) .* (H + beta .* H_der(n, beta .* r)) .* P_der .* cos(m .* ph);
                Eph = - 1./ (beta .* r) .* (H + beta .* H_der(n, beta .* r)) .* m .* P(m, 1)./sin(theta) .* sin(m .* ph);
            end
        end
        
        Erho_s(i, :, :) = Eth .* cos(theta) + Er .* sin(theta);
        Ephi_s(i, :, :) = Eph;
        Ez_s(i, :, :) = Er .* cos(theta) - Eth .* sin(theta);
                
        Hrho_s = - Ephi_s./Z0;
        Hphi_s = Erho_s./Z0;
        
end

%% Normalization of free space

Q_fs = Pow(1:1:M, Erho_s, Ephi_s, Hrho_s, Hphi_s, rho, drho, dphi);

%% Coupling between waveguide and free space

for k = 1:1:N(end)
    for l = 1:1:M
        X_i = (squeeze(Ep_rho(k, :, :)) .* squeeze(Hphi_s(l, :, :)) - squeeze(Ep_phi(k, :, :)) .* squeeze(Hrho_s(k, :, :))) .* rho .* drho .* dphi;
        X(k, l) = sum(sum(X_i));
    end
end

Ip = eye(length(M), length(M));
Ir = eye(length(N), length(N));

F_ = 2 * inv(Q_wg + X * inv(Q_fs) * X.');

Spp = inv(Q_fs) * X.' * F_ * X - Ip;
Spr = inv(Q_fs) * X.' * F_ * Q_wg;
% Srp = 2 * eye(length(Nr), length(Np)) - F_ * X;
%Srp = F_ * X;
Srp = Spr.';
Srr = F_ * Q_wg - Ir;

S = [Spp Spr; Srp Srr];

%% 

r = 1000 .* lamb;
theta = -pi-eps:dphi:pi-eps;
ph = eps:dphi:2*pi;
[theta, ph] = meshgrid(theta, ph);

Eth = zeros(M, size(theta, 1), size(theta, 2));
Eph = zeros(M, size(theta, 1), size(theta, 2));
Er = zeros(M, size(theta, 1), size(theta, 2));

M = 10;
Z0 = 120 * pi;

for i = 1:1:M
    
        S = load('Sper.mat');
        Sper = S.Sper;
        
        mode = Sper(i).mode;
        pol = Sper(i).pol;
        n = Sper(i).n;
        m = Sper(i).m;
    
        H = besselh(n, 2, beta .* r);
        P = legendre(n, cos(theta));
        P_der = legendre_der_(n, m, P, cos(theta));
        
        if mode == "TE"
            if pol == "even"
                Er(i, :, :) = zeros(size(theta)) + eps;
                Eth(i, :, :) = H .* m .* P(m, 1) ./ sin(theta) .* sin(m .* ph);
                Eph(i, :, :) = sin(theta) .* H .* m .* P_der .* cos(m .* ph);
            else
                Er(i, :, :) = zeros(size(theta)) + eps;
                Eth(i, :, :) = H .* m .* P(m, 1) ./ sin(theta) .* cos(m .* ph);
                Eph(i, :, :) = sin(theta) .* H .* m .* P_der .* sin(m .* ph);
            end
        end
        
        
        if mode == "TM"
            if pol == "even"
                Er(i, :, :) = n .* (n + 1) .* H ./ (beta .* r) .* P(m, 1) .* sin(m .* ph);
                Eth(i, :, :) = - sin(theta) .* 1./ (beta .* r) .* (H + beta .* H_der(n, beta .* r)) .* P_der .* sin(m .* ph);
                Eph(i, :, :) = - 1./ (beta .* r) .* (H + beta .* H_der(n, beta .* r)) .* m .* P(m, 1)./sin(theta) .* cos(m .* ph);
            else
                Er(i, :, :) = n .* (n + 1) .* H ./ (beta .* r) .* P(m, 1) .* cos(m .* ph);
                Eth(i, :, :) = - sin(theta) .* 1./ (beta .* r) .* (H + beta .* H_der(n, beta .* r)) .* P_der .* cos(m .* ph);
                Eph(i, :, :) = - 1./ (beta .* r) .* (H + beta .* H_der(n, beta .* r)) .* m .* P(m, 1)./sin(theta) .* sin(m .* ph);
            end
        end
        
end



ar = ones(N, 1);
bp = Spr * ar;
Eth_t = zeros(size(theta));
Eph_t = zeros(size(theta));
Er_t = zeros(size(theta));

for n = 1:1:M
    Eth_t = Eth_t + squeeze(Eth(n, :, :)) .* bp(n);
    Eph_t = Eph_t + squeeze(Eph(n, :, :)) .* bp(n);
    Er_t = Er_t + squeeze(Er(n, :, :)) .* bp(n);
end


E = sqrt(abs(Er_t).^2 + abs(Eth_t).^2 + abs(Eph_t).^2);
figure;
plot(theta(1, :)*180/pi, db(abs(E(1, :))/max(abs(E(1, :)))), 'LineWidth', 2);
hold on;
plot(theta(91, :)*180/pi, db(abs(E(91, :))/max(abs(E(91, :)))), 'LineWidth', 2);
grid on;


