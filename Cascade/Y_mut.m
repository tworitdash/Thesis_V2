function [Ymn_mut] = Y_mut(m, n, x, y, a, b, dx, dy, k, omega, mu, epsilon)

    C = 1j .* 8 * a ./ (b .* omega .* mu);
    TermY = b - y;
    TermX = fmn(x, m, n, a, k);
    
    Ph = exp(-1j .* k .* sqrt(x.^2 + y.^2))./(sqrt(x.^2 + y.^2));
    
    if m == n
        Ymn_i = C .* TermY .* TermX .* Ph;
    else
        Ymn_i = C .* TermY .* TermX .* Ph;
    end
    Ymn_mut = sum(sum(Ymn_i)) .* dx .* dy;
    
%     if (m == 5) && (n == 5)
%         hold on;
%        figure;
%         surface(x, y, abs(Ymn_i)); shading flat;
%         colorbar;
%         hold on;
%     end
end