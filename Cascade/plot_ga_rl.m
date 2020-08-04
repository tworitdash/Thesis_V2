
F_feko = linspace(4e9, 8e9, 36);


F = 4e9:0.5e9:21e9;

data = read(rfdata.data,'../../Thesis/3wg/gsm/RL_Ga.s2p');
s_params = extract(data,'S_PARAMETERS');

RL_feko = db(abs(squeeze(s_params(1, 1, :)) + squeeze(s_params(1, 1, :))).^2)./2;

figure(2);
hold on;

plot(F_feko * 1e-9, RL_feko, 'LineWidth', 2); grid on;
hold on;

