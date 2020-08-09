function [Ep_rho, Ep_phi] = E_ap(N, rho, phi, F, r, z, epsilon, mu)

        c0 = 3e8;
        Str = load('Xmn.mat');
        Xmn = Str.Xmn;
        
        for i = N
            
            disp(i);
            mode = Xmn(i).mode;
            m = Xmn(i).m;
            xmn = Xmn(i).xmn;
            
            beta_rho = xmn./r;
            beta = (2 .* pi .* F)./c0;
            
            omega = 2 .* pi .* F;
            
            beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
            
            if m == 0
                deltan = 1;
                
            else
                deltan = 0;
            end
            
            if mode == "TE"
                A = sqrt(2 ./ pi) ./ (besselj(m, xmn) .* sqrt(xmn.^2 - m.^2) .* sqrt(1 + deltan));
%                 A = 1;
                Ph = exp(-1j .* beta_z .* z);
                Ep_rho(i, :, :) = A.^2 .* m ./ rho .* besselj(m, beta_rho .* rho) .* sin(m .* phi) .* Ph;
                Ep_phi(i, :, :) = A.^2 .* beta_rho .* besselj_der(m, beta_rho .* rho) .* cos(m .* phi) .* Ph;
            elseif mode == "TM"
                B = sqrt(2 ./ pi) ./ (xmn .* besselj(m + 1, xmn) .* sqrt(1 + deltan));
%                 B = 1;
                Ph = exp(-1j .* beta_z .* z);
                Ep_rho(i, :, :) = B.^2 .* beta_rho .* besselj_der(m, beta_rho .* rho) .* sin(m .* phi) .* Ph;
                Ep_phi(i, :, :) = B.^2 .* m ./ rho .* besselj(m, beta_rho .* rho) .* cos(m .* phi) .* Ph;
            end

        end
        
end