function [STT, STR, SRT, SRR, N_] = GSM_N_SameAzimuth(R, L, er, mur, F, k)
c0 = 3e8;

n = length(R);

Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;
XMN = horzcat(str.xmn);

for i = 1:n
    f =  fc(R(i), er(i), mur(i));
%     N_i  =  find(f < F);
    xmn_i = 2.*pi.*R(i)./c0 .* f(end);
    idx = find(XMN <= xmn_i);
    N_(i) = length(idx);
end

if (N_(end) < 50)
    %N = [20 round(20 * (R(end)/R(1))^2)];
    N = N_ + k;
else
    N = N_;
end


% N = round(linspace(5, 25, n));

% N = [20 20 20] ; % Number of modes

J = length(R) - 1; % Number of Junctions
%% Frequency independent inner cross product 

for j = 1:J
    x_til = zeros(N(j), N(j + 1));
    x_til(:, :) = Inner_p2_SameAzimuth(1:1:N(j), 1:1:N(j + 1), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j));
    X_til(j).x_til = x_til;
end

%% Frequency loop to find the GSM of the entire structure

for k = 1:length(F)
    
    disp('Frequency Iteration: ');
    disp(k);

if n == 2
    
    [STT_, STR_, SRT_, SRR_] = GSM_V2_SameAzimuth(1:1:N(1), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), X_til(1).x_til);

else

[S33, S34, S43, S44] = GSM_V2_SameAzimuth(1:1:N(1), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), X_til(1).x_til);
[S11, S12, S21, S22] = GSM_V2_SameAzimuth(1:1:N(2), 1:1:N(3), F(k), R(3), R(2), er(3), mur(3), er(2), mur(2), X_til(2).x_til);
Sl = SL(R(2), F(k), 1:1:N(2), L(2));

  
[STT_, STR_, SRT_, SRR_] = cascade_3(1:1:N(2), S11, S12, S21, S22, S33, S34, S43, S44, Sl);

% Use the for loop in case of more than 3 junctions (J > 3)

for j = 3:J

    % recursion 
    
    [S11, S12, S21, S22] = GSM_V2_SameAzimuth(1:1:N(j), 1:1:N(j + 1), F(k), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j), X_til(j).x_til);
    S33 = STT_; S34 = STR_; S43 = SRT_; S44 = SRR_;
    Sl = SL(R(j), F(k), 1:1:N(j), L(j));
    
    [STT_, STR_, SRT_, SRR_] = cascade_3(1:1:N(j), S11, S12, S21, S22, S33, S34, S43, S44, Sl);
    
end

end

slr = SL(R(1), F(k), 1:1:N(1), L(1));
slt = SL(R(end), F(k), 1:1:N(end), L(end));

STT(k, :, :) = slt * STT_ * slt'; 
STR(k, :, :) = slt * STR_ * slr; 
SRT(k, :, :) = slr * SRT_ * slt; 
SRR(k, :, :) = slr * SRR_ * slr;


end

end