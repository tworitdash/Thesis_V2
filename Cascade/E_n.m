function [Erho_n, Ephi_n] = E_n(Nr, rho, phi, F, r, z, epsilon, mu, drho, dphi)
        c0 = 3e8;
        Str = load('Xmn.mat');
        Xmn = Str.Xmn;
        
        
        
        for i = 1:length(Nr)
            disp(i);
            mode = Xmn(i).mode;
            m = Xmn(i).m;
            xmn = Xmn(i).xmn;
            
            beta_rho = xmn./r;
            beta = (2 .* pi .* F)./c0;
            
            omega = 2 .* pi .* F;
            
            if mode == "TE"
%                 
                    A = 1;
                    C1 = 1; % C1 and D1 are for the rho component
                    D1 = 0; 

                    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
                    
                    [Erho, Ephi, ~] = E_TE(epsilon, m, rho, phi, beta_rho, z, beta);
                    [Hrho, Hphi, ~] = H_TE(epsilon, m, rho, phi, beta_rho, z, beta, omega, mu);

%                     Erho = -A .* m./(epsilon .* rho) .* besselj(m, beta_rho .* rho) .* (-C1 .* sin(m .* phi)...
%                         + D1 .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
%                     Ephi = A .* beta_rho./epsilon .* besselj_der(m, beta_rho .* rho) .* (C1 .* cos(m .* phi)...
%                         + D1 .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
%                     Ez = ones(size(rho)) .* 1e-5;
                    
                    
                    
                    
                    
%                     Hrho = -A .* (beta_rho .* beta_z) ./ (omega .* mu .* epsilon) .* besselj_der(m, beta_rho .* rho) .* (C1 .* cos(m .* phi)...
%                         + D1 .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
%                     Hphi = -A .* (m .* beta_z ./ (omega .* mu .* epsilon .* rho)) .* besselj(m, beta_rho .* rho) .* (-C1 .* sin(m .* phi)...
%                         + D1 .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
                    
                    Power = P(Erho, Ephi, Hrho, Hphi, rho, drho, dphi);
                    
                    Erho_n(i, :, :) = Erho; % ./ sqrt(Power);
                    Ephi_n(i, :, :) = Ephi; % ./ sqrt(Power);
                    
%                     E_i(i, :, :) = sqrt(Erho.^2 + Ephi.^2 + Ez.^2);
                
            elseif mode == "TM"
                
                
                   B = 1;
                   C = 1;
                   D = 0;
    
                    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
                    
                    [Erho, Ephi, ~] = E_TM(epsilon, m, rho, phi, beta_rho, z, beta, omega, mu);
                    [Hrho, Hphi, ~] = H_TM(mu, m, rho, phi, beta_rho, z, beta);
    
%                     Erho = -B .* (beta_rho .* beta_z ./ (omega .* mu .* epsilon)) .* besselj_der(m, beta_rho .* rho).* (C .* cos(m .* phi)...
%         + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
%                     Ephi = -B .* (m .* beta_z ./ (omega .* mu .* epsilon .* rho)) .* besselj(m, beta_rho .* rho) .* (- C .* sin(m .* phi)...
%         + D .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
% %                     Ez(i, :, :) = -1j .* B .* (beta_rho.^2 ./ (omega .* mu .* epsilon)) .* besselj(m, beta_rho .* rho) .*  (C .* cos(m .* phi)...
% %         + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
% 
%                     Hrho = -A .* m./(mu .* rho) .* besselj(m, beta_rho .* rho) .* (-C1 .* sin(m .* phi)...
%                         + D1 .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
%                     Hphi = A .* beta_rho./mu .* besselj_der(m, beta_rho .* rho) .* (C1 .* cos(m .* phi)...
%                         + D1 .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
                    
                    Power = P(Erho, Ephi, Hrho, Hphi, rho, drho, dphi);
                    
                    Erho_n(i, :, :) = Erho; % ./ sqrt(Power);
                    Ephi_n(i, :, :) = Ephi; % ./ sqrt(Power);
%                    
            end
            
            
        end
        
%         E = E_i;
    
end