function [Power] = P(N, Erho, Ephi, Hrho, Hphi, rho, drho, dphi)

   for i = 1:length(N)
        P_i = (squeeze(Erho(i, :, :)) .* (squeeze(Hphi(i, :, :))) - (squeeze(Hrho(i, :, :))) .* squeeze(Ephi(i, :, :))) .* rho .* drho .* dphi;
    
        Power(i, i) = sum(sum(P_i));
   end

end