function [RLRR_TE11] = GSM_N_opt_allvar_V2_freq(R_cone, Len, F, k)

c0 = 3e8;
lamb = c0./5e9;

% R1 = R_vec(1);
% Rend = R_vec(2);

% R = linspace(R1, Rend, round(Len/(lamb/10)));

% R = linspace(R1, Rend, E);
lamb_opt_freq = c0./F(end);

num = round(Len./(lamb_opt_freq./10));

n_R = length(R_cone);

N_axis = round(num./(n_R-1));

for p = 1:n_R-1
    R_(:, p) = linspace(R_cone(p), R_cone(p+1), N_axis);
end

R = reshape(R_, 1, size(R_, 1) .* size(R_, 2));

n = length(R);

l1 = lamb/4;
L = [l1 ones(1, n - 1) * Len/n];

er = ones(n);
mur = ones(n);

for i = 1:n
    f =  fc(R(i), er(i), mur(i));
    N_i  =  find(f < F(end));
    N_(i) = length(N_i);
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
    x_til(:, :) = Inner_p2(1:1:N(j), 1:1:N(j + 1), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j));
    X_til(j).x_til = x_til;
end

%% Frequency loop to find the GSM of the entire structure
RLRR_TE11 = size(zeros(1, length(F)));

parfor k = 1:length(F)
    
    disp('Frequency Iteration: ');
    disp(k);

if n == 2
    
    [STT_, STR_, SRT_, SRR_] = GSM_V2(1:1:N(1), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), X_til(1).x_til);

else

[S33, S34, S43, S44] = GSM_V2(1:1:N(1), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), X_til(1).x_til);
[S11, S12, S21, S22] = GSM_V2(1:1:N(2), 1:1:N(3), F(k), R(3), R(2), er(3), mur(3), er(2), mur(2), X_til(2).x_til);
Sl = SL(R(2), F(k), 1:1:N(2), L(2));

  
[STT_, STR_, SRT_, SRR_] = cascade_3(1:1:N(2), S11, S12, S21, S22, S33, S34, S43, S44, Sl);

% Use the for loop in case of more than 3 junctions (J > 3)

for j = 3:J

    % recursion 
    
    [S11, S12, S21, S22] = GSM_V2(1:1:N(j), 1:1:N(j + 1), F(k), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j), X_til(j).x_til);
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

RLRR_TE11(k) = db(sum(sum(abs(SRR(k, :, :)).^2)))./2; % Return loss at waveguide R
end


end