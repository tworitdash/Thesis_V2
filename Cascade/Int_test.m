dth = pi/1800;

theta = 1e-2/2:dth:2*pi;

k = 1:100;

for i = 1:length(k)

Y = sin(theta) .* cos(k(i) .* abs(cos(theta))) * dth;

Y_int(i) = sum(Y);
end

plot(k, Y_int*1e16);
% 
% Y2 = 2 .* pi .* StruveH0(k)
%% 

rho = linspace(1e-7, a, 1000);
drho = rho(2) - rho(1);

k = 10;
l = 30;

Int_i = besselj(0, abs(k .* rho)) .* exp(-1j .* l .* rho);
Int = sum(Int_i) .* drho

figure;
plot(rho, real(Int_i));
hold on;
plot(rho, imag(Int_i));

%% 
m = 3;
n = 5;
rho = 1:1:1000;
z1 = n * pi/a .* rho;
z2 = m * pi/a .* rho;
J = StruveH0(z1) - StruveH0(z2);

plot(z, J);



