c0 = 3e8;

F = 14e9;
lamb = c0/F;

N = 13;

R = 2.3e-2;

z = 0.005; % + lamb/4;

er = 1;
mur = 1;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er .* er0;
mu = mur .* mu0;


[rho, phi] = meshgrid(eps:R/100:R, eps:pi/180:2*pi+eps);


[Ep_rho, Ep_phi] = E_ap(1:1:N, rho, phi, F, R, z, epsilon, mu);


E_rho = zeros(size(rho));
E_phi = zeros(size(rho));


for i = 1:1:N
    
%     figure;
    
    EM_rho = conj(squeeze(Ep_rho(i, :, :))); %  * exp(1j .* pi/2);
    EM_phi = conj(squeeze(Ep_phi(i, :, :))); % *  exp(1j .* pi/2);
        
    E_rho = E_rho + EM_rho;
    E_phi = E_phi + EM_phi;
   
end

E_abs_full = sqrt(abs(E_rho).^2 + abs(E_phi).^2);

% figure;
% 
% surface(rho .* cos(phi), rho .* sin(phi), db(abs(E_rho)./max(max(abs(E_rho))))); shading flat;
% 
% % max(max(abs(E_rho)))
% 
% colormap('jet');
% 
% figure;
% 
% surface(rho .* cos(phi), rho .* sin(phi), db(abs(E_phi)./max(max(abs(E_phi))))); shading flat;
% 
% colormap('jet');



figure;

surface(rho .* cos(phi), rho .* sin(phi), db(abs(E_rho))); shading flat;

% max(max(abs(E_rho)))

colormap('jet');

figure;

surface(rho .* cos(phi), rho .* sin(phi), db(abs(E_phi))); shading flat;

colormap('jet');


% figure;
% 
% surface(rho .* cos(phi), rho .* sin(phi), (real(E_rho))); shading flat;
% title('Re E_{\rho}')
% % max(max(abs(E_rho)))
% 
% colormap('jet');
% 
% figure;
% 
% surface(rho .* cos(phi), rho .* sin(phi), (real(E_phi))); shading flat;
% title('Re E_{\phi}')
% colormap('jet');
% 
% 
% figure;
% 
% surface(rho .* cos(phi), rho .* sin(phi), (imag(E_rho))); shading flat;
% title('Im E_{\rho}')
% % max(max(abs(E_rho)))
% 
% colormap('jet');
% 
% figure;
% 
% surface(rho .* cos(phi), rho .* sin(phi), (imag(E_phi))); shading flat;
% title('Im E_{\phi}')
% colormap('jet');

Plot_NF_feko_V4;
