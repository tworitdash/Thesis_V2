function [Ymn_mut] = ymn(m, n, a, b, x, y, k0, mur, er, dx, dy)
    c0 = 3e8;
    v = c0./sqrt(er);
    
    omega = k0 .* v;
    
    mu0 = 1.25663706e-6;  % Free Space Permeability
    mu = mu0 .* mur;
    
   
    
   C = 1j .* 8 .* a ./ (b .* omega .* mu);
%    C = 1j;
   
   ph = exp(-1j .* k0 .* sqrt(x.^2 + y.^2))./(sqrt(x.^2 + y.^2));
   
   [x_term] = fmn_rec(m, n, a, x, k0);
   [y_term] = b - y;
   
   Ymn_i = C .* y_term .* x_term .* ph .* dx .* dy;
   Ymn_mut = sum(sum(Ymn_i));
end