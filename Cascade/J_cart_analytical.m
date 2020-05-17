function [J_ft_x, J_ft_y, J_ft_z] = J_cart_analytical(E0, k0, theta_obs, phi_obs, a, b, M, Dm, Gamma)
    

    kx = k0 .* sin(theta_obs) .* cos(phi_obs);
    ky = k0 .* sin(theta_obs) .* sin(phi_obs);
    kz = k0 .* cos(theta_obs);
    
    C0 = b .* sin(ky .* b./2)./(ky .* b./2);
    C1 = (2 * 1 .* pi .* a .* (1j).^(1 - 1) .* cos(kx .* a ./ 2))./((1*pi).^2 - (kx .* a).^2);
    
    CmDm = zeros(size(theta_obs));
    
    for i = 1:M
        Cm = (2 * M(i) .* pi .* a .* (1j).^(M(i) - 1) .* cos(kx .* a ./ 2))./((M(i)*pi).^2 - (kx .* a).^2); 
        CmDm = CmDm + Dm(i) .* Cm;
    end
    
    J_ft_x = 1e-7;
    J_ft_y = E0 .* (1 + Gamma) .* C0 .* (C1 + CmDm);
    J_ft_z = 1e-7;
   
 
end