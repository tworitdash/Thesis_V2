%% Zeros of Bessel's function (For TM) and zeros of the derivative of the Bessel's function (For TE) Calculation
clear;

m = 1:1:100;
% m = 1:1:100;
% m = 1;

for j = 1:length(m)

disp(j);

xm = linspace(0.1, 1000, 10000);

Jm = @(z) besselj(m(j), z);  % Bessel's function

dz = 1e-5;

Jm_der = @(z) (besselj(m(j), z + dz) - besselj(m(j), z - dz))./(2 * dz); % Derivative of the Bessel's function

ym_TE = Jm_der(xm); % Values of the Bessel's function at the positions defined by xm
ym_TM = Jm(xm);

chsign_TE = find(diff(sign(ym_TE)));  % Detection of the sign changes 
chsign_TM = find(diff(sign(ym_TM)));

J_TE = size(chsign_TE, 2); % number of sign changes
J_TM = size(chsign_TM, 2);

for i = 1:J_TE
    
    xmn_TE((j - 1) * J_TE + i).xmn = fzero(Jm_der, xm(chsign_TE(i)));  % finding the roots near the points where Jm changes sign
    xmn_TE((j - 1) * J_TE + i).m = m(j);
    xmn_TE((j - 1) * J_TE + i).n = i;
    xmn_TE((j - 1) * J_TE + i).mode = "TE";
    disp(i);
end

for k = 1:size(chsign_TM, 2)
    xmn_TM((j - 1) * J_TM + k).xmn = fzero(Jm, xm(chsign_TM(k)));
    xmn_TM((j - 1) * J_TM + k).m = m(j);
    xmn_TM((j - 1) * J_TM + k).n = k;
    xmn_TM((j - 1) * J_TM + k).mode = "TM";
    disp(k);

end

end

Xmn = [xmn_TE xmn_TM];



[x,idx]=sort([Xmn.xmn]);

Xmn = Xmn(idx);

[x_TE, idx_TE] = sort([xmn_TE.n]);
xmn_TE = xmn_TE(idx_TE);

[x_TE, idx_TE] = sort([xmn_TE.m]);
xmn_TE = xmn_TE(idx_TE);

% 
% [x,idx]=sort([Xmn.n]);
% Xmn = Xmn(idx);
% 
% [x,idx]=sort([Xmn.m]);
% Xmn = Xmn(idx);




% save('Xmn', 'Xmn');
save('Xmn_azimuthal_inc_TE', 'xmn_TE');

