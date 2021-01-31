function [Power] = Pow_2modes(Erho, Ephi, Hrho, Hphi, rho, drho, dphi)

 
        P_i = (Erho .* Hphi - Hrho .* Ephi) .* rho .* drho .* dphi;
    
        Power = sum(sum(P_i));

end