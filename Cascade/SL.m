function [SL] = SL(rr, F, Nr, L)

    c0 = 3e8;
    
    
    Str = load('Xmn.mat');
    Xmn = Str.Xmn;
   
    SL = zeros(length(Nr), length(Nr));
    
    for i = 1:length(Nr)
        beta = (2 * pi * F) ./ c0;
        beta_rho = (Xmn(i).xmn)./rr;
        beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
        SL(i, i) = exp(-1j .* beta_z .* L);
    end
end