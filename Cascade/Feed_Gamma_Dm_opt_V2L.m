function [ModeNumberAper, Transmission_sum]  = Feed_Gamma_Dm_opt_V2L(R, L, F, k, er, mur)
 

c0 = 3e8;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability


omega = 2 * pi * F;

lamb = c0./F;

k0 = 2.*pi./lamb;


% Length = 5e-2;
% Length = 0.5 .* lamb;


% drho = R(end)/100;
% dphi = pi/180;

% [rho, ph] = meshgrid(eps:drho:R(end), eps:dphi:2*pi+eps);

epsilon = er .* er0;
mu = mur .* mu0;

% L = ones(1, n) .* Length/n; % length of each waveguide section

% L(1) = 2.5 * L(1);
% L(end) = 1.25e-2;
% L(1) = L(1) + lamb./4;

% [E_aperture_rho, E_aperture_phi] = Near_field_fun_2(er, mur, R, F, L, rho, ph, 20, drho, dphi);


[STT, STR, SRT, SRR, N] = GSM_N_SameAzimuth(R, L, er, mur, F, k);

Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;
%% 
ModeNumberAper = N(end);


ap = zeros(N(end), 1);
ar = zeros(N(1), 1);
ar(1) = 1;

if N(1) == 1
    STR_req = STR(1, 1:N(end));
else
    STR_req = squeeze(STR(1, 1:N(end), 1:N(1)));
end

bp =  STR_req * ar;  %+ squeeze(STT(1, 1:N(end), 1:N(end))) * ap;
Transmission_sum = bp;
% 
% %%
% 
% for i = 1:ModeNumberAper
%     HigherModes = i+1:1:i+HM;
%     
%     for d = 1:length(HigherModes)
%     
%         xmn(d) = str(HigherModes(d)).xmn;
%         beta_rho(d) = xmn(d)./R(end);
%         M(d) = str(HigherModes(d)).m;
%         beta_z(d) = -1j .* sqrt(-(k0.^2 - (beta_rho(d)).^2));
%         YTE(d) = beta_z(d)./(omega .* mu(end));
%         ZTE(d) = 1./YTE(d);
%     
%     end
%     
%  for l = 1:length(HigherModes)
%     for p = 1:length(HigherModes)
%         disp('iteration inside');
%         disp(p);
%         disp('iteration outside');
%         disp(l);
%         
%         if l == p
%             Ymut(l, p) = Yin_Circular(HigherModes(l), HigherModes(p), k0, R(end), er(end), mur(end), timesk0) + YTE(p);
%         else
%             Ymut(l, p) = Yin_Circular(HigherModes(l), HigherModes(p), k0, R(end), er(end), mur(end), timesk0);
%         end
%     end
% end
% 
% for e = 1:length(HigherModes)
%     Y_rhs(e) = Yin_Circular(i, HigherModes(e), k0, R(end), er(end), mur(end), 100);
% end
% 
% Dm(i, :) = Ymut\(-Y_rhs.');
% 
% Yii = Yin_Circular(i, i, k0, R(end), er(end), mur(end), 100);
% 
% 
% xmnii = str(i).xmn;
% beta_rhoii = xmnii./R(end);
% 
% beta_zii = -1j .* sqrt(-(k0.^2 - beta_rhoii.^2));
% YTEii = beta_zii./(omega .* mu(end));
% ZTEii = 1./YTEii;
% 
% yii(i) = Yii./YTEii;
% 
% yap(i) = (Yii + Dm(i, :) * Y_rhs.')./YTEii;
% 
% Gamma(i) = (1 - yap(i))./(1 + yap(i));
% end
% 
% %% Far fields
% % [theta, phi] = meshgrid(eps:dtheta:pi+eps, eps:dphi:2*pi+eps);
% 
% % Eth = zeros(size(theta));
% % Eph = zeros(size(theta));
% % 
% % for o = 1:ModeNumberAper
% %     
% %     HigherModes = o+1:1:o+20;
% %     
% %     [Eth_o, Eph_o] = FF_apertureFSCir(o, length(HigherModes)+1, [1, Dm(o, :)], Gamma(o), theta, phi, F, er(end), mur(end), R(end));
% %     Eth = Eth + Eth_o .* Transmission_sum(o);
% %     Eph = Eph + Eph_o .* Transmission_sum(o);
% % 
% % end
% 
% % E_FF = sqrt(abs(Eth).^2 + abs(Eph).^2);
% 
% %% plots
% % figure;
% % 
% % plot(theta(1, :).*180/pi, db(E_FF(1, :)), 'LineWidth', 2);
% % hold on;
% % plot(theta(91, :).*180/pi, db(E_FF(91, :)), 'LineWidth', 2);
% 
% % grid on;
% % ylim([-50 30]);
% 
% %% Directivity
% 
% % zeta = 120 .* pi;
% % 
% % U = abs(E_FF).^2./(2 .* zeta);
% % 
% % P_rad_i = U(:, 91:end) .* sin(theta(:, 91:end)) .* dtheta .* dphi;
% % 
% % P_rad = sum(sum(P_rad_i));
% % 
% % D = 4 .* pi .* U ./ P_rad;
% 
% % 
% % figure;
% % plot(theta(1, :), db(D(1, :))/2, 'LineWidth', 2);
% % hold on;
% % plot(theta(91, :), db(D(91, :))/2, 'LineWidth', 2);
% % 
% % grid on;
end
