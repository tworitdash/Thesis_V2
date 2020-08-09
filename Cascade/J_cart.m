function [J_ft_x, J_ft_y, J_ft_z] = J_cart(J_cx, J_cy, J_cz, dx, dy, k0, x, y, z, theta_obs, phi_obs)
    
    
     for i = 1:size(theta_obs, 1)
         for j = 1:size(theta_obs, 2)
             
            kx = k0 .* sin(theta_obs(i, j)) .* cos(phi_obs(i, j));
            ky = k0 .* sin(theta_obs(i, j)) .* sin(phi_obs(i, j));
            kz = k0 .* cos(theta_obs(i, j));

            J_i_x = J_cx .* exp(1j .* (kx .* x + ky .* y));
            J_i_y = J_cy .* exp(1j .* (kx .* x + ky .* y));
            J_i_z = J_cz .* exp(1j .* (kx .* x + ky .* y));

            J_ft_x(i, j) = sum(sum(J_i_x)) * dx * dy;
            J_ft_y(i, j) = sum(sum(J_i_y)) * dx * dy;
            J_ft_z(i, j) = sum(sum(J_i_z)) * dx * dy;
         end
     end
 
end