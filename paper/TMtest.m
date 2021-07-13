c0 = 3e8;
er = 1;
mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

% f = linspace(10.15e9, 14e9, 41);
f = linspace(4.99654e9, 20e9, 41);
k0 = 2 * pi .* f./c0;
omega = 2 .* pi .* f;


Str = load('Xmn_azimuthal_inc_TM.mat');

str = Str.xmn_TM;

% N = linspace(1, 10, 10);
N = 1;

R = 1.8e-2;

for i = 1:length(N)
    for k = 1:length(f)
        xmn(i, k) = str(N(i)).xmn;
        beta_rho(i, k) = xmn(i)./R;
        M(i, k) = str(N(i)).m;
        beta_z(i, k) = -1j .* sqrt(-(k0(k).^2 - (beta_rho(i, k)).^2));
        ZTM(i, k) = beta_z(i, k)./(omega(k) .* epsilon);
        YTM(i, k) = 1./ZTM(i, k);
       
    end
     figure(1000); hold on; plot(f*1e-9, abs(ZTM(i, :)), 'LineWidth', 2);
end
%%
c0 = 3e8;
er = 1;
mur = 1;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;
% f = linspace(4.99654e9, 20e9, 41);

k0 = 2 * pi .* f./c0;
omega = 2 .* pi .* f;


Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;

% N = linspace(1, 10, 10);
N = 1;

R = 1.8e-2;

for i = 1:length(N)
    for k = 1:length(f)
        xmn(i, k) = str(N(i)).xmn;
        beta_rho(i, k) = xmn(i)./R;
        M(i, k) = str(N(i)).m;
        beta_z(i, k) = -1j .* sqrt(-(k0(k).^2 - (beta_rho(i, k)).^2));
        YTE(i, k) = beta_z(i, k)./(omega(k) .* mu);
        
       
    end
     figure(1001); hold on; plot(f*1e-9, abs(YTE(i, :)), 'LineWidth', 2);
end