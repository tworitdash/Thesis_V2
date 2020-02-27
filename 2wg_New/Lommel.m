function [X] = Lommel(a, b, beta_rhoNu, beta_rhoMu, Nu, Mu)


    if (beta_rhoNu == beta_rhoMu) && (Nu == Mu)
        Xa = 1./2 .* a^2 .* (besselj(Nu, beta_rhoNu .* a).^2 - ...
besselj(Nu -1, beta_rhoNu .* a) .*  besselj(Nu + 1, beta_rhoNu .* a)); 
        Xb = 1./2 .* b^2 .* (besselj(Nu, beta_rhoNu .* b).^2 - ...
besselj(Nu -1, beta_rhoNu .* b) .*  besselj(Nu + 1, beta_rhoNu .* b)); 
        X = Xb - Xa;
    elseif (beta_rhoNu ~= beta_rhoMu) && (Nu == Mu)
        Xa = Zeta(a, beta_rhoNu, beta_rhoMu, Nu, Mu);
        Xb = Zeta(b, beta_rhoNu, beta_rhoMu, Nu, Mu);
    
        X = Xb - Xa;
    else
        X = 0;
    end

end