function [RLRR_TE11, SRR, STT, STR, SRT] = GSM_N_opt_allvar_V2_freq_new(R_cone, Len, F, k)

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

r_len = length(R_(1:end-1, 1));

for q = 1:n_R-1
    for r = 1:r_len
       R((q - 1) .* r_len + r) = R_(r, q);
    end
end

R = [R  R_(end, end)];

% R = reshape(R_, 1, size(R_, 1) .* size(R_, 2));

n = length(R);

l1 = lamb/4;
L = [l1 ones(1, n - 1) * Len/(n - 1)];

for i = 1:n
    L_axis(i) = sum(L(1:i));
end
L_axis = [0 L_axis];
figure(1); hold on;
plot(L_axis, [R(1) R], 'o', 'LineWidth', 2, 'Color', [0, 0.0780, 0.1840]);
hold on;
plot(L_axis, [-R(1) -R], 'o', 'LineWidth', 2, 'Color', [0, 0.0780, 0.1840] );
title(['Optimum Horn, Aprture Radius = ', num2str(R(end)./lamb), ' \lambda ', 'Horn Length = ', num2str(sum(L)./lamb), ' \lambda'], 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Horn Length cut', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Horn Radius cut', 'FontSize', 12, 'FontWeight', 'bold');
grid on;



er = ones(n);
mur = ones(n);

for i = 1:n
    f =  fc(R(i), er(i), mur(i));
    N_i  =  find(f < F(end));
    N_(i) = length(N_i);
end

if (N_(end) < 100)
    %N = [20 round(20 * (R(end)/R(1))^2)];
    N = N_ + k;
else
    N = N_;
end

% N = round(linspace(5, 25, n));

% N = [20 20 20] ; % Number of modes

J = length(R) - 1; % Number of Junctions
%% Frequency independent inner cross product 

parfor j = 1:J
    disp('Junction');
    disp(j);
    disp('of')
    disp(J);
   
    x_til = zeros(N(j), N(j + 1));
    x_til(:, :) = Inner_prod(1:1:N(j), 1:1:N(j + 1), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j));
    X_til(j).x_til = x_til;
end

%% Frequency loop to find the GSM of the entire structure
RLRR_TE11 = zeros(1, length(F));

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

STT(k, :, :) = slt * STT_ * slt; 
STR(k, :, :) = slt * STR_ * slr; 
SRT(k, :, :) = slr * SRT_ * slt; 
SRR(k, :, :) = slr * SRR_ * slr;
SRR_1 = squeeze(SRR(k, :, :));

f_base =  fc(R(1), er(1), mur(1));
Num_modes_prop  =  find(f_base < F(k));

% RLRR_TE11(k) = db(sum(sum(abs(SRR(k, 1:Num_modes_prop(end), 1:Num_modes_prop(end))).^2)))./2; % Return loss at waveguide R
RLRR_TE11(k) = db(abs(sum(SRR_1(1, 1:Num_modes_prop(end)))).^2)./2; % Return loss at waveguide R
% RLRR_TE11(k) = db(abs(sum(SRR(k, 1, 1:Num_modes_prop(end))))).^2./2; % Return loss at waveguide R
end



end