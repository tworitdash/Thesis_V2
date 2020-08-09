clear;
a = 8.3e-2;
b = 1e-2;
M = 2 * (1:15) + 1;


c0 = 3e8;
F = 1.6e9;
omega = 2* pi * F;
k = omega./c0;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

mu = mu0;
epsilon = er0;


for i = 1:length(M)
    for j = 1:length(M)
        m = M(i);
        n = M(j);
        if m == n
            ep = 1e-6;
            
            Y1(i, j) = Y_mut2(m, n, ep, ep, k, omega, mu, epsilon);
            
            dx = a/100;
            dy = b/100;
            
            [x_, y_] = meshgrid(eps:dx:ep, eps:dy:ep);
            
            Y12(i, j) = Y_mut(m, n, x_, y_, ep, ep, dx, dy, k, omega, mu, epsilon);
            
            x = ep:dx:a;
            y = ep:dy:b;
            
            [x, y] = meshgrid(x, y);
            
            Y2(i, j) = Y_mut(m, n, x, y, a, b, dx, dy, k, omega, mu, epsilon);
            
            Y(i, j) = Y1(i, j) + Y2(i, j);

        else
            x = eps:dx:a;
            y = eps:dy:b;

            [x, y] = meshgrid(x, y);
            
            Y(i, j) = Y_mut(m, n, x, y, a, b, dx, dy, k, omega, mu, epsilon);
            
        end
      
%         figure;
       
%         surface(x, y, abs()); shading flat;
    end
end
% figure;
% for i = 1:length(M)
%     for j = 1:length(M)
%         
%         if i == j
%             figure(2)
%             hold on;
%             plot(x, squeeze(a(i, j, :)));
%             hold on;
%         else
%             figure(1);
%             hold on;
%             plot(x, squeeze(a(i, j, :)), 'color', [0, i/20, 1 - i/100]);
%             hold on;
%         end
%        
%         
%         
%     end
% end



