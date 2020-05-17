function [Ymn] = Y_mut2(m, n, a, b, k, omega, mu, epsilon)

    r = sqrt(a^2 + b^2);
    
    rho = linspace(eps, r, 10000);
    
%     tau = 25;
    
    drho = rho(2) - rho(1);

    C = 1j .* 8 * a ./ (b .* omega .* mu);
    
    if m == n
        C1 = 1./(4 .* pi^2 .* a .* n);
        zeta1 = (k.^2 + (n .* pi./a).^2);
        zeta2 = (k.^2 - (n .* pi./a).^2);
        Konst = C .* C1 .* 2 .* b;
        Int_rho = (zeta1 .* StruveH0(n * pi * rho./a) + zeta2 .* pi .* n .* besselj(0, abs(n * pi * rho./a))) .* exp(-1j .* k .* rho);
        Ymn = Konst .* sum(Int_rho) .* drho;
    else
        C1 = -(1j).^(m + n)./(2 * pi^2 * a * (n^2 - m^2));
        zeta1 = (k.^2 - (m .* pi ./a).^2) .* n ;
        zeta2 = (k.^2 - (n .* pi ./a).^2) .* m ;
        Konst = C .* C1 .* 2 .* pi .* b;
        Int_rho = (zeta1 .* StruveH0(m .* pi .* rho./a) - zeta2 .* StruveH0(n * pi .* rho./a)) .* exp(-1j .* k .* rho);
        Ymn =  Konst .* sum(Int_rho) .* drho;
    end
    
%     if (n == 17) && (m == 3)
%         figure(3);
%         hold on;
%         plot(rho*1e3, abs(Int_rho), 'Linewidth', 2);
%         hold on;
%         grid on;
%     end
end