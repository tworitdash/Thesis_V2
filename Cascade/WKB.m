%% parameters for the horn

c0 = 3e8;

er = 1; mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

F = 5e9;
omega = 2 * pi * F;
beta = omega / c0;

r0 = 2e-2;

r1 = 4e-2;

z = linspace(eps, 4.5e-2, 2);

slope = (r1 - r0)/z(end);

r = r0 + z .* slope;




% plot(z, r); % taper

%% Modes

Str = load('Xmn.mat');
Xmn = Str.Xmn;



for k = 1:length(z)
    
%     [fc_] = fc(r(k), er, mur);
% 
%      N = find(fc_ < F);
       N = 1;
     [rho, phi] = meshgrid(eps:r(k)/100:r(k),  eps:pi/180:2*pi-eps);
     
     E_z_sum = zeros(size(rho));
     E_rho_sum = zeros(size(rho));
     E_phi_sum = zeros(size(rho));

for i = 1:1:N(end)
    
    
    
    xmn = Xmn(i).xmn;
    m = Xmn(i).m;
    n = Xmn(i).n;
    mode = Xmn(i).mode;
    
    
    beta_rho = xmn./r(k);
    
    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    
    
    Psi_1(i, k) = acos(xmn./(beta .* r(1)));
    Psi_end(i, k) = acos(xmn./(beta .* r(k)));
    
    Psi = [Psi_1(i, k) Psi_end(i, k)];
    
    if mode == "TE"

%         K = beta_z ./ (omega .* mu .* epsilon^2);
       
        
        Z_i = 2 * pi * F * mu./ beta_z;
        Y_i = 1./Z_i;
        K = 1j.*beta_z./beta_rho.^2 .* Z_i;
        
    elseif mode == "TM"

%         K = beta_z ./ (omega .* mu.^2 .* epsilon);
        Z_i = beta_z ./ (2 * pi * F .* epsilon);
        Y_i = 1./Z_i;
        
        
        K = 1j.*beta_z./beta_rho.^2;
        
    end
    
    if mode == "TM"
        h = (2./sqrt(beta_z)) .* exp(-1j .* xmn .*((tan(Psi(2)) - Psi(2))) - (tan(Psi(1)) - Psi(1)));
        Pot =  besselj(m, beta_rho .* rho) .* cos(m .* phi) .* h;
        
        E_z_i(i, :, :) = Pot .* h; 
        E_z_sum = E_z_sum + squeeze(E_z_i(i, :, :));
        
        E_rho_i(i, :, :) = -K .* beta_rho .* besselj_der(m, beta_rho .* rho) .* cos(m .* phi) .* h;
        E_rho_sum = E_rho_sum + squeeze(E_rho_i(i, :, :));
        
        E_phi_i(i, :, :) = K .* m./rho .* besselj(m, beta_rho .* rho) .* sin(m .* phi) .* h;
        E_phi_sum = E_phi_sum + squeeze(E_phi_i(i, :, :));
    
    else 
        h = (2./sqrt(beta_z)) .* exp(-1j .* xmn .*((tan(Psi(2)) - Psi(2))) - (tan(Psi(1)) - Psi(1)));
        Pot =  besselj(m, beta_rho .* rho) .* cos(m .* phi); 
        
        E_z_i(i, :, :) = zeros(size(rho, 1), size(rho, 2));
        E_z_sum = E_z_sum + squeeze(E_z_i(i, :, :));
        
        E_rho_i(i, :, :) = K .* m./rho .* besselj(m, beta_rho .* rho) .* sin(m .* phi) .* h;
        E_rho_sum = E_rho_sum + squeeze(E_rho_i(i, :, :));
        
        E_phi_i(i, :, :) = K .* beta_rho .* besselj_der(m, beta_rho .* rho) .* cos(m .* phi) .* h;
        E_phi_sum = E_phi_sum + squeeze(E_phi_i(i, :, :));
        
        
    end
    
end

% E_z_k(k, :, :) = squeeze(sum(E_z_i, 1));
E_z_k(k, :, :) = E_z_sum;
E_rho_k(k, :, :) = E_rho_sum;
E_phi_k(k, :, :) = E_phi_sum;
end


%% Plot

x = rho.*cos(phi);
y = rho.*sin(phi);
E_z = squeeze(E_z_k(end, :, :));
E_rho = squeeze(E_rho_k(end, :, :));
E_phi = squeeze(E_phi_k(end, :, :));

figure;
surface(rho.*cos(phi), rho.*sin(phi), db(abs(E_z)./(max(max(abs(E_z)))))); shading flat;
% surface(rho.*cos(phi), rho.*sin(phi), db(abs(E_z))); shading flat;
colormap('jet');

figure;
surface(rho.*cos(phi), rho.*sin(phi), db(abs(E_rho)./(max(max(abs(E_rho)))))); shading flat;
% surface(rho.*cos(phi), rho.*sin(phi), db(abs(E_z))); shading flat;
colormap('jet');

figure;
surface(rho.*cos(phi), rho.*sin(phi), db(abs(E_phi)./(max(max(abs(E_phi)))))); shading flat;
% surface(rho.*cos(phi), rho.*sin(phi), db(abs(E_z))); shading flat;
colormap('jet');