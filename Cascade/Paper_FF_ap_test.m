clear;
close all;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

R = 3e-2;
er = 1; mur = 1;

f = fc(R, er, mur);

F = f(1)+f(1)/10;

N = 1;
Dm = 1;
Gamma = 0;
dth = pi/180;
dph = pi/180;
[theta, phi] = meshgrid(eps:dth:pi/2, eps:dph:2*pi);

[Eth, Eph, Eco, Exp, CO, XP] = FF_ap_test(N, Dm, Gamma, theta, phi, F, er, mur, R);

Eres = sqrt(abs(Eth).^2 + abs(Eph).^2);

figure; plot(theta(1, :)*180/pi, db(abs(Eres(1, :))));