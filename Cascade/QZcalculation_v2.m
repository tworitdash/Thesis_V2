%% Normalization matrix with different modes on the waveguide:
function [Q, Z, Y, K] = QZcalculation_v2(N, F, r, er, mur, rho_, phi_, z, drho, dphi)
c0 = 3e8;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur;


% Q = zeros(size(F, 2), N(end), N(end));
Q = zeros(length(N), length(N));
Z = zeros(length(N), length(N));
Y = zeros(length(N), length(N));


xmn_ = zeros(length(N));
% beta_z_ = zeros(N);

K = zeros(length(N), length(N));


%Y = zeros(m(end), m(end));

omega = 2 * pi * F;

lamb = c0./F; % wavelength
 
beta = 2 * pi ./ lamb;



for i = 1:length(N)
    
    disp(i);
    
       Str = load('Xmn.mat');
       Xmn = Str.Xmn;
      
      
      xmn_(i) = Xmn(N(i)).xmn;
      m = Xmn(N(i)).m;
      
      mode = Xmn(N(i)).mode;
      
if m == 0
%% Numerical Q     

[Erho, Ephi, Ez, Hrho, Hphi, Hz, beta_z] = E_and_H_v2(rho_, phi_, er, mur, z, r, N(i), F);
   

    
    Poyn = (Erho .* Hphi - Hrho .* Ephi) .* rho_ * drho .* dphi;
    Qij = sum(sum(Poyn));
    
    
    
    if mode == "TE"

        K(i, i) = beta_z ./ (omega .* mu .* epsilon^2);
        
    elseif mode == "TM"

        K(i, i) = beta_z ./ (omega .* mu.^2 .* epsilon);
        
    end
%     
%     if mode == "TE"
%             Norm = (epsilon * pi/2 .* (xmn_(i).^2 - m.^2) .* (besselj(m, xmn_(i))).^2).^(-1);
%         elseif mode == "TM"
%             Norm = (epsilon .* pi/2 .* xmn_(i).^2 .* (besselj_der(m, xmn_(i))).^2).^(-1);
%     end
    Norm = 1;
    
    Q(i, i) = Qij .* Norm;




else
    
 %% Analytical Q
 
%       
      fc = xmn_(i) ./ (2 * pi * r * sqrt(mu .* epsilon));
      disp(fc);
      
      beta_rho = xmn_(i)./r;

      beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    
      if mode == "TE"

        K(i, i) = beta_z ./ (omega .* mu .* epsilon^2);
        
      elseif mode == "TM"

        K(i, i) = beta_z ./ (omega .* mu.^2 .* epsilon);
        
      end
%       
%                 if mode == "TE"
%                     Norm = (epsilon * pi/2 .* (xmn_(i).^2 - m.^2) .* (besselj(m, xmn_(i))).^2).^(-1);
%                 elseif mode == "TM"
%                     Norm = (epsilon .* pi/2 .* xmn_(i).^2 .* (besselj_der(m, xmn_(i))).^2).^(-1);
%                 end
        Norm = 1;        
               
      Const = K(i, i) .* beta_rho.^2 ./ 4 .* Norm;  

        A = Lommel(0, r, beta_rho, beta_rho, m - 1, m - 1);

        C = Lommel(0, r, beta_rho, beta_rho, m + 1, m + 1);
        
%         if m == 0
%             Isin = 0;
%             Icos = 2 * pi;
%             Qij = (Isin + Icos) * beta_rho.^4 .* Lommel(0, r, beta_rho, beta_rho, 1, 1);
%         else
            Isin = intphisin(0, 2*pi, m, m);
            Icos = intphicos(0, 2*pi, m, m);
            Qij = Const .* ((Isin + Icos) .* (A + C));
%         end
    
      
    
       Q(i, i) = Qij;
end 
%      
%% Impedance (Z) and Admittance (Y) 
      
    
    if mode == "TE"
        Z_i = 2 * pi * F * mu./ beta_z;
        Y_i = 1./Z_i;
    elseif mode == "TM"
        Z_i = beta_z ./ (2 * pi * F .* epsilon);
        Y_i = 1./Z_i;
    end
    

    
    Z(i, i) = Z_i;
    Y(i, i) = Y_i;
end


end