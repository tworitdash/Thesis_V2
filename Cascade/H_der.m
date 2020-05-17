function [H_d] = H_der(n, z)
        
        H_d = 1./2 .* (besselh(n - 1, 2, z) - besselh(n + 1, 2, z));

end