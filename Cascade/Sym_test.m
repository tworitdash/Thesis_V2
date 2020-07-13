k0 = 2*pi*5e9/3e8;
th_ = -pi/2-eps:pi/180:pi/2-eps;
beta_rhop = 1.8;
arg = k0 .* sin(th_);

R = 4e-2;
for i = 1:length(th_)
    I0(i) = Lommel(0, R, beta_rhop, arg(i), 0, 0);
    I2(i) = Lommel(0, R, beta_rhop, arg(i), 2, 2);
end

figure;
plot(th_*180/pi, I0);
figure;
plot(th_*180/pi, I2);