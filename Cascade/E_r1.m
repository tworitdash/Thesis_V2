function [Erho, Ephi, Ez] = E_r(Nr, rho, phi, F, r, z, epsilon, mu)
        c0 = 3e8;
        Str = load('Xmn.mat');
        Xmn = Str.Xmn;
        
        E_i = zeros(length(Nr), size(rho, 1), size(rho, 2));
        er = 1;
        mur = 1;
        
        z_ = 0;
        R = 2.3e-2;
        
        drho = R./100; 
        dphi = pi/180;
        
        [Q, Z, Y, K] = QZcalculation_v2(Nr, F, R, er, mur, rho, phi, z_, drho, dphi);
        
        for i = 1:length(Nr)
            disp(i);
            mode = Xmn(i).mode;
            m = Xmn(i).m;
            xmn = Xmn(i).xmn;
            
            beta_rho = xmn./r;
            beta = (2 .* pi .* F)./c0;
            
            omega = 2 .* pi .* F;
            
            if m == 0
                    deltan = 1;
            else
                    deltan = 0;
            end
            
            if mode == "TE"
                
                    
                    A = K(i, i) .* Z(i, i) .* sqrt(abs((1 + deltan).*(pi/2 .* (xmn.^2 - m.^2) .* (besselj(m, xmn)).^2).^(-1)));
%                      A = sqrt(abs((pi/2 .* (xmn.^2 - m.^2) .* (besselj(m, xmn)).^2).^(-1)));
%                     A = 1;
   
                    C1 = 1; % C1 and D1 are for the rho component
                    D1 = 0; 

                    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));

                    Erho(i, :, :) = - A .* m./(rho) .* besselj(m, beta_rho .* rho) .* (-C1 .* sin(m .* phi)...
                        + D1 .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
                    Ephi(i, :, :) = A .* beta_rho .* besselj_der(m, beta_rho .* rho) .* (C1 .* cos(m .* phi)...
                        + D1 .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
                    Ez(i, :, :) = ones(size(rho)) .* 1e-5;
                    
%                     E_i(i, :, :) = sqrt(Erho.^2 + Ephi.^2 + Ez.^2);
                
            elseif mode == "TM"
                
                   B = K(i, i) .* Z(i, i) .* sqrt(abs((1 + deltan).*(pi/2 .* (xmn).^2 .* (besselj_der(m, xmn)).^2).^(-1)));
%                    B = sqrt(abs((pi/2 .* (xmn).^2 .* (besselj_der(m, xmn)).^2).^(-1)));
%                    B = 1;
   
                   C = 1;
                   D = 0;
    
                    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    
                    Erho(i, :, :) = -B .* (beta_rho) .* besselj_der(m, beta_rho .* rho).* (C .* cos(m .* phi)...
        + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
                    Ephi(i, :, :) = -B .* (m ./ (rho)) .* besselj(m, beta_rho .* rho) .* (- C .* sin(m .* phi)...
        + D .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
                    Ez(i, :, :) = -1j .* B .* (beta_rho.^2 ./ (omega .* mu .* epsilon)) .* besselj(m, beta_rho .* rho) .*  (C .* cos(m .* phi)...
        + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
                    
%                     E_i(i, :, :) = sqrt(Erho.^2 + Ephi.^2 + Ez.^2);
            end
            
            
        end
        
%         E = E_i;
    
end