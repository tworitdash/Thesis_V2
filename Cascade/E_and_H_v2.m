function [Erho, Ephi, Ez, Hrho, Hphi, Hz, beta_z] = E_and_H_v2(rho_, phi_, er, mur, z, r, n, f)

%% Waveguide parameter


c0 = 3e8;
er0 = 8.85418782e-12;
mu0 = 1.25663706e-6;

epsilon = er * er0;
mu = mu0 * mur;
%% Zeros of Bessel's functions and their derivatives:
    
    Str = load('Xmn.mat');
    
    Xmn = Str.Xmn;
    m = Xmn(n).m;
    mode = Xmn(n).mode;
    xmn = Xmn(n).xmn;

    

    beta_rho = xmn./r;  % wave number along the rho direction (\beta_{\rho})

%% Cut of frequency

fc = xmn ./ (2 * pi * r * sqrt(mu .* epsilon));
disp(fc);

%% Source Parameters
% f = fc + eps;

omega = 2 * pi * f;

lamb = c0./f; % wavelength
 
beta = 2 * pi ./ lamb; % wave number

beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));

%% Wave equations

if mode == "TE"
    
    [Erho, Ephi, Ez] = E_TE(epsilon, m, rho_, phi_, beta_rho, z, beta);
    [Hrho, Hphi, Hz] = H_TE(epsilon, m, rho_, phi_, beta_rho, z, beta, omega, mu);
   
elseif mode == "TM"
    [Erho, Ephi, Ez] = E_TM(epsilon, m, rho_, phi_, beta_rho, z, beta, omega, mu);
    [Hrho, Hphi, Hz] = H_TM(mu, m, rho_, phi_, beta_rho, z, beta);
end



end