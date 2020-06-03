clear;

c0 = 3e8;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

F = 5e9;

lamb = c0./F;

rr = 2e-2; % Base radius
% rt = 4e-2; % Top redius
rt = rr .* 4;

n = 150; % number of transitions

Length = 5e-2;

R = linspace(rr, rt, n); % radius vector

drho = R(end)/100;
dphi = pi/180;

[rho, ph] = meshgrid(eps:drho:R(end), eps:dphi:2*pi+eps);

n = length(R);

er = ones(1, n); % Relative Permittivity of each WG section
mur = ones(1, n); % Relative Permeability of each WG sectio

epsilon = er .* er0;
mu = mur .* mu0;

L = ones(1, n) .* Length/n; % length of each waveguide section

L(1) = 2.5 * L(1);
L(end) = 1.25e-2;


% [E_aperture_rho, E_aperture_phi] = Near_field_fun_2(er, mur, R, F, L, rho, ph, 20, drho, dphi);


[STT, STR, SRT, SRR, N] = GSM_N_SameAzimuth(R, L, er, mur, F, 0);

Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;

ModeNumberAper = N(end);

HigherModes = N(end)+1:1:N(end)+20;

for i = 1:length(HigherModes)
    
    xmn(i) = str(HigherModes(i)).xmn;
    beta_rho(i) = xmn(i)./R;
    M(i) = str(HigherModes(i)).m;
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


    
    





