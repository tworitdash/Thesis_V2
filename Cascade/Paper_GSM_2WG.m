%% GSM 2 waveguide junction
clear;
close all;

c0 = 3e8;
f_l = 5e9;
lamb = c0/f_l;


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability


R = [1.8e-2 1.8e-2*3];

L = ones(1, length(R)) .* lamb/4;

er = ones(1, length(R));
mur = ones(1, length(R));


F = linspace(5e9, 15e9, 20);


for i = 1:length(R)
    f =  fc(R(i), er(i), mur(i));
    N_i  =  find(f < F(end));
    N(i) = length(N_i);
end

x_til = Inner_p2(1:1:N(2), 1:1:N(2), R(2), R(1), er(2), mur(2), er(1), mur(1));


parfor k = 1:length(F)

    [STT_, STR_, SRT_, SRR_] ...
        = GSM_V2(1:1:N(2), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), x_til);
    
    slr = SL(R(1), F(k), 1:1:N(end), L(1));
    slt = SL(R(end), F(k), 1:1:N(end), L(end));

    STT(k, :, :) = slt * STT_ * slt'; 
    STR(k, :, :) = slt * STR_ * slr; 
    SRT(k, :, :) = slr * SRT_ * slt; 
    SRR(k, :, :) = slr * SRR_ * slr;

end


%% Plots
data5_feko = read(rfdata.data,'../../Thesis/Paper_Thesis/Paper_2WG.s20p');
s_params_5_feko = extract(data5_feko,'S_PARAMETERS');


figure;
hold on;
plot(F * 1e-9, db(abs(squeeze(SRR(:, 1, 1)))), 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(s_params_5_feko(1, 1, :)))), 'LineWidth', 2); grid on;

xlim([5 15]);


% %% Surface Plots
% 
% figure; 
% imagesc(1:10, F * 1e-9, db(abs(squeeze(SRT(:, 1, 1:10)))));shading flat; colormap('jet')
% 
% figure; 
% imagesc(1:10, F * 1e-9, db(abs(squeeze(s_params_5_feko(1, 11:20, :)))));shading flat; colormap('jet')