function [Erho, Ephi, Ez, Hrho, Hphi, Hz] = E_TE_FEKOstyle(mu, m, rho, phi, beta_rho, z, beta)
    
    c0 = 3e8;

    A = 1;
   
    C1 = 1; % C1 and D1 are for the rho component
    D1 = 0; 
%     
%     C2 = 0; % C2 and D2 are for the phi component
%     D2 = 1;
    omega = beta .* c0;
    
    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
   
    ZTE = (omega .* mu ./ beta_z);
    
    Erho = -1j .* A .* m .* beta_z ./(beta_rho.^2 .* rho) .* ZTE .* besselj(m, beta_rho .* rho) .* (-C1 .* sin(m .* phi)...
        + D1 .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
    Ephi = 1j .* A .* beta_z ./ beta_rho  .* ZTE .* besselj_der(m, beta_rho .* rho) .* (C1 .* cos(m .* phi)...
        + D1 .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
    Ez = zeros(size(rho)) + eps;
    
    Hrho = -Ephi/ZTE;
    Hphi = Erho/ZTE;
    
    Hz = A .* besselj(m, beta_rho .* rho) .* cos(m .* phi);
end