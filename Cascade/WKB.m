%% parameters for the horn

c0 = 3e8;

er = 1; mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

F = 14e9;
omega = 2 * pi * F;
beta = omega / c0;

r0 = 2e-2;

r1 = 4e-2;

z = linspace(eps, 5e-2, 1000);

slope = (r1 - r0)/z(end);

r = r0 + z .* slope;




% plot(z, r); % taper

%% Modes

Str = load('Xmn.mat');
Xmn = Str.Xmn;



for k = 1:length(z)
    
    [fc_] = fc(r(k), er, mur);

     N = find(fc_ < F);

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

        K = beta_z ./ (omega .* mu .* epsilon^2);
        
        Z_i = 2 * pi * F * mu./ beta_z;
        Y_i = 1./Z_i;
        
    elseif mode == "TM"

        K = beta_z ./ (omega .* mu.^2 .* epsilon);
        Z_i = beta_z ./ (2 * pi * F .* epsilon);
        Y_i = 1./Z_i;
        
    end
    
    
    [rho, phi] = meshgrid(eps:r(k)/100:r(k),  eps:pi/180:2*pi-eps);
    
    if mode == "TM"
    
        E_z_i(i, :, :) = sqrt(K) .* sqrt(Z_i) .* (beta_rho.^2./beta_z) .* besselj(m, beta_rho .* rho) .* cos(m .* phi) .* (2./sqrt(beta_z)) .* ...
            exp(-1j .* xmn .*((tan(Psi(2)) - Psi(2))) - (tan(Psi(1)) - Psi(1))); 
    
    else 
        
        E_z_i(i, :, :) = zeros(size(rho, 1), size(rho, 2));
    end
    
end

E_z_k(k, :, :) = squeeze(sum(E_z_i, 1));

end


%% Plot

x = rho.*cos(phi);
y = rho.*sin(phi);
E_z = squeeze(E_z_k(end, :, :));


figure;
surface(rho.*cos(phi), rho.*sin(phi), db(abs(E_z)./(max(abs(E_z))))); shading flat;
colormap('jet');
