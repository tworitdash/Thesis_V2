clear;
c0 = 3e8; % speed of light
R = 2e-2; % Radius of the waveguide

er = 1;
mur = 1;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

% F = 5e9;
% lamb = c0/F;

x = linspace(0.6, 1, 41);
lamb = 2.*R./x;

f = c0./lamb;

% r = lamb .* x/2;

% lamb = c0/F;


L = 100;

N = 2:1:20;

Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;

for b = 1:length(f)
    
%     R = r(b);
    
      F = f(b);
      omega = 2 * pi * F;

      k0 = omega./c0;

for i = 1:length(N)
    
    xmn(i) = str(N(i)).xmn;
    beta_rho(i) = xmn(i)./R;
    M(i) = str(N(i)).m;
    beta_z(i) = -1j .* sqrt(-(k0.^2 - (beta_rho(i)).^2));
    YTE(i) = beta_z(i)./(omega .* mu);
    ZTE(i) = 1./YTE(i);
    
end

for l = 1:length(N)
    for p = 1:length(N)
        disp('iteration inside');
        disp(p);
        disp('iteration outside');
        disp(l);
        
        if l == p
            Ymut(l, p) = Yin_Circular(N(l), N(p), k0, R, er, mur, L) + YTE(p);
        else
            Ymut(l, p) = Yin_Circular(N(l), N(p), k0, R, er, mur, L);
        end
    end
end

for l = 1:length(N)
    Y_rhs(l) = Yin_Circular(1, N(l), k0, R, er, mur, L);
end

Dm(b, :) = Ymut\(-Y_rhs.');

Y11 = Yin_Circular(1, 1, k0, R, er, mur, L);


xmn11 = str(1).xmn;
beta_rho11 = xmn11./R;

beta_z11 = -1j .* sqrt(-(k0.^2 - beta_rho11.^2));
YTE11 = beta_z11./(omega .* mu);
ZTE11 = 1./YTE11;

y11(b) = Y11./YTE11;

yap(b) = (Y11 + Dm(b, :) * Y_rhs.')./YTE11;

Gamma(b) = (1 - yap(b))./(1 + yap(b));
Gamma_TE11(b) = (1 - y11(b))./(1 + y11(b));

end

figure(1);
hold on;
plot(x, real(yap), 'LineWidth', 2); grid on; hold on; plot(x, imag(yap), 'LineWidth', 2);
hold on;
plot(x, real(y11), '-.', 'LineWidth', 1); grid on; hold on; plot(x, imag(y11), '-.', 'LineWidth', 1);
% 
figure(3);
hold on;
plot(x, db(abs(Gamma)), 'LineWidth', 2); grid on;
figure(4);
plot(x, angle(Gamma)*180/pi, 'LineWidth', 2); grid on;

% figure;
% surf(x, N, db(abs(Dm.'))); shading flat;
% xlabel('2 r/\lambda', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('n in TE_{1, n} at the aperture', 'FontSize', 12, 'FontWeight', 'bold');
% title('Higher order mode excitations at the aperture (D_m) in dB');
% colorbar;
% colormap('jet');

figure(3);
hold on;
plot(x, db(abs(Gamma_TE11)), 'LineWidth', 2); grid on;
figure(4);
hold on;
plot(x, angle(Gamma_TE11)*180/pi, 'LineWidth', 2); grid on;