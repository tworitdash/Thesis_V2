function [Spp, Spr, Srp, Srr] = GSM(Nr, Np, F, rp, rr, erp, murp, err, murr, X_)

Spp = zeros(length(Np), length(Np));
Spr = zeros(length(Np), length(Nr));
Srp = zeros(length(Nr), length(Np));
Srr = zeros(length(Nr), length(Nr));

%% Frequency independent Modular inner cross product between the two wavegudies

X_til = zeros(length(Nr), length(Np));

 for p = 1:length(Np)
     for r = 1:length(Nr)
           X_til(r, p) = X_(r, p);
     end
 end
    
    

%% Wavwguide p

drho = rp/100;
dphi = pi/180;

[rho_, phi_] = meshgrid(eps:drho:rp, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zp = 0; 

[Qp, Zp, Yp, Kp] = QZcalculation_v2(Np, F, rp, erp, murp, rho_, phi_, zp, drho, dphi);


%% Wavwguide r

drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

[Qr, Zr, Yr, Kr] = QZcalculation_v2(Nr, F, rr, err, murr, rhor_, phir_, zr, drho, dphi);


Ip = eye(length(Np), length(Np));
Ir = eye(length(Nr), length(Nr));

X = sqrt(Kr * Zr) * X_til * sqrt(Yp * Kp); % modular inner cross product. Takes the dimension of Np \times Nr

F_ = 2 * inv(Qr + X * inv(Qp) * X.');

Spp = inv(Qp) * X.' * F_ * X - Ip;
Spr = inv(Qp) * X.' * F_ * Qr;
% Srp = 2 * eye(length(Nr), length(Np)) - F_ * X;
%Srp = F_ * X;
Srp = Spr.';
Srr = F_ * Qr - Ir;

% M = Qp + X.' * inv(Qr) * X;
% N = X.' * inv(Qr) * X - Qp;
% 
% Srr = -2 * inv(Qr) * X * inv(M) * X.' + Ir;
% Srp = inv(Qr) * X * (Ip - inv(M) * N);
% Spr = 2 * inv(M) * X.';
% Spp = inv(M) * N;

end


