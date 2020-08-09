function [Erho, Ephi, Ez, Hrho, Hphi, Hz] = E_TM_FEKOstyle(epsilon, m, rho, phi, beta_rho, z, beta, omega, mu)
    c0 = 3e8;
    B = 1;
   
    C = 1;
    D = 0;
   
    omega = beta .* c0;
    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    ZTM = beta_z./(omega .* epsilon);
    
    Erho = -1j .*B * beta_z./beta_rho .* besselj_der(m, beta_rho .* rho).* (C .* cos(m .* phi)...
        + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
    Ephi = -1j .* B .* (m .* beta_z ./ (beta_rho.^2 .* rho)) .* besselj(m, beta_rho .* rho) .* (- C .* sin(m .* phi)...
        + D .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
    Ez = -1j .* B .* besselj(m, beta_rho .* rho) .*  (C .* cos(m .* phi)...
        + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z) ;
    
    Hrho = -Ephi./ZTM;
    Hphi = Erho./ZTM;
    Hz = zeros(size(rho)) + eps;
end