
delGx = @(x, y) (-x .* exp(-1j .* sqrt(x.^2 + y.^2))) .* (1j .* sqrt(x.^2 + y.^2) - 1) ./ (x.^2 + y.^2).^(3./2);

delGy = @(x, y) (-y .* exp(-1j .* sqrt(x.^2 + y.^2))) .* (1j .* sqrt(x.^2 + y.^2) - 1) ./ (x.^2 + y.^2).^(3./2);

ep1 = 2e-6;
ep2 = 1e-5;

Int1 = integral(@(t) delGy(ep1, t), -ep2, ep2);
Int2 = integral(@(t) delGx(t, ep2), ep1, -ep1);
Int3 = integral(@(t) delGy(-ep1, t), ep2, -ep2);
Int4 = integral(@(t) delGx(t, -ep2), -ep1, ep1);

Int = Int1 + Int2 + Int3 + Int4

%% 

x = linspace(-1e-5, 1e-5, 100);
y = linspace(-1e-5, 1e-5, 100);

[x, y] = meshgrid(x, y);

delGx = (x .* exp(-1j .* sqrt(x.^2 + y.^2))) .* (1j .* sqrt(x.^2 + y.^2) - 1) ./ (x.^2 + y.^2).^(3./2);

delGy = (y .* exp(-1j .* sqrt(x.^2 + y.^2))) .* (1j .* sqrt(x.^2 + y.^2) - 1) ./ (x.^2 + y.^2).^(3./2);

figure;
%quiver(x, y, abs(delGx), abs(delGy), 'LineWidth', 5);

surface(x, y, abs(sqrt(delGx.^2 + delGy.^2))); shading flat; colormap('jet');


%% 

x = linspace(-1e-5, 1e-5, 100);
y = linspace(-1e-5, 1e-5, 100);

[x, y] = meshgrid(x, y);

Laplacian =  (exp(-1j .* sqrt(x.^2 + y.^2))) .* (1j.^2 .* (x.^2 + y.^2) + 1j .* sqrt(x.^2 + y.^2) + 1) ./ (x.^2 + y.^2).^(3./2);

sum(sum(Laplacian))
